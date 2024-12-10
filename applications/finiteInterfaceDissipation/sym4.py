import sympy as sp
import numpy as np

def replace_powers_with_intpow(expr):
    """
    Recursively replace power expressions in the form x**n with intpow<n>(x) using SymPy.
    
    Args:
        expr (sp.Basic): The sympy expression to process.

    Returns:
        sp.Basic: The modified expression as a sympy object.
    """
    def transform_pow(expr):
        if isinstance(expr, sp.Pow):  # Check if it's a power operation
            base, exp = expr.args
            # Only replace integer exponents and ensure it's not a division case
            if exp.is_Integer and True:#exp >= 0:
                return sp.Function(f"dealii::Utilities::fixed_power<{int(exp)}>")(transform_pow(base))
                #return sp.Function(f"intpow<{int(exp)}>")(transform_pow(base))
        
        # Recursively apply transformation to arguments
        if expr.args:
            return expr.func(*[transform_pow(arg) for arg in expr.args])
        
        return expr

    # Transform the entire expression
    return transform_pow(expr)

def int_to_float(expr):
    return expr.xreplace({n: sp.Float(n) for n in expr.atoms(sp.Rational)})

def wrap_decimal_numbers(expr):
    """
    Wraps decimal numbers in constV() in a sympy expression tree.
    
    Args:
        expr (sp.Basic): The sympy expression to process.

    Returns:
        sp.Basic: The modified expression with decimal numbers wrapped.
    """
    def transform_numbers(expr):
        def transform_if_num(expr):
            if expr.is_Number:  # Check for floating-point literals
                return sp.Function("constV")(expr)
            else:
                return transform_numbers(expr)

        if isinstance(expr, sp.Add):
            return sum(transform_if_num(arg) for arg in expr.args)
        # Recursively apply transformation to arguments
        if expr.args:
            return expr.func(*[transform_numbers(arg) for arg in expr.args])
        return expr
    
    return transform_numbers(expr)

def partial_derivatives_to_cpp(expr, variables):
    """
    Takes a symbolic expression, computes its partial derivatives,
    and formats them for C++ code generation.

    Args:
        expr (sp.Expr): The symbolic expression.
        variables (list): A list of sympy symbols with respect to which the partial derivatives are taken.

    Returns:
        dict: A dictionary where keys are the variables and values are the formatted C++ strings.
    """
    # Compute partial derivatives
    derivatives = {var: sp.diff(expr, var) for var in variables}

    # Apply transformations
    derivative_strings = {}
    for var, derivative in derivatives.items():
        # Replace powers with intpow and wrap decimals in constV()
        transformed_expr = replace_powers_with_intpow(derivative)
        transformed_expr = int_to_float(transformed_expr)
        transformed_expr = wrap_decimal_numbers(transformed_expr)
        transformed_expr = int_to_float(transformed_expr)
        
        # Convert the expression to a string suitable for C++
        cpp_formatted_str = str(transformed_expr)
        derivative_strings[str(var)] = cpp_formatted_str

    return derivative_strings

def L_at_temp(L_arr, T):
    return [L_arr[0]+L_arr[1]*T, L_arr[2]+L_arr[3]*T]

if __name__ == "__main__":
    # 1:Ga
    # 2:In
    # 3:Pb
    pairname = {'12':'Ga-In', '23':'In-Pb', '31':'Pb-Ga'}


    GaIn0a, GaIn0b, GaIn1a, GaIn1b,\
    InPb0a, InPb0b, InPb1a, InPb1b,\
    PbGa0a, PbGa0b, PbGa1a, PbGa1b = sp.symbols('GaIn0a, GaIn0b, GaIn1a, GaIn1b,\
                                                 InPb0a, InPb0b, InPb1a, InPb1b,\
                                                 PbGa0a, PbGa0b, PbGa1a, PbGa1b')


    L_data = {'Ga-In' : np.array([ GaIn0a, GaIn0b, GaIn1a, GaIn1b]),
              'In-Pb' : np.array([ InPb0a, InPb0b, InPb1a, InPb1b]),
              'Pb-Ga' : np.array([ PbGa0a, PbGa0b, PbGa1a, PbGa1b])
        }
    
    #L_data = {'Ga-In' : np.array([ 5.14148219e+03, -7.89869349e-01,  3.14491130e+03, -9.80783643e+00]),
    #          'In-Pb' : np.array([ 3.775e+03, -1.285e+00,  1.830e+02,  3.810e-01                    ]),
    #          'Pb-Ga' : np.array([ 2.28142140e+04, -8.99292499e+00, -4.90642270e+03, 5.96724182e+00 ])
    #    }


    ref_H = {'Ga': [0.0, 20000, 20000],
              'In(A6)': [9000, 0.0, 4644.0],
              'alpha': [20000, 38.0, 967.0],
              'Pb(A1)': [20000, 42.0, 0.0]}
    ref_S = {'Ga': [0.0, 0.0, 0.0],
              'In(A6)': [0.0, 0.0, 0.233],
              'alpha': [0.0, 0.002, 1.183],
              'Pb(A1)': [0.0, 0.012, 0.0]}
    L_terms = {'Ga': [0.0, 0.0, 0.0, 0.0],
                'In(A6)': [-2859, 0, 3341, 0],
                'alpha': [3246, 0, 939, 0],
                'Pb(A1)': [5069, -2.624, 456, 0.875]
                }

    # Define symbols
    T = sp.symbols('T')
    R = sp.symbols('R')
    x_A, x_B, x_C = sp.symbols('x_A.val, x_B.val, x_C.val')
    # Define variables for partial differentiation
    variables = [x_A, x_B, x_C]

    x1 = x_A
    x2 = x_B
    x3 = x_C
    # 
    x3 = 1.-x2-x1

    H = ref_H['alpha']
    S = ref_S['alpha']

    #solid
    g_ref = x1*(H[0]-T*S[0]) + x2*(H[1]-T*S[1]) + x3*(H[2]-T*S[2])
    #liquid
    #g_ref = 0.0

    g_ent = R*T*(x1*sp.log(x1)+x2*sp.log(x2)+x3*sp.log(x3))
    
    G_excess12 = (L_at_temp(L_data[pairname['12']], T)[0]+(x1-x2)*L_at_temp(L_data[pairname['12']], T)[1])
    G_excess23 = (L_at_temp(L_data[pairname['23']], T)[0]+(x2-x3)*L_at_temp(L_data[pairname['23']], T)[1])
    G_excess31 = (L_at_temp(L_data[pairname['31']], T)[0]+(x3-x1)*L_at_temp(L_data[pairname['31']], T)[1])

    # Example expression
    W12 = 4.*x1*x2/(1.-(x1-x2)*(x1-x2))
    W23 = 4.*x2*x3/(1.-(x2-x3)*(x2-x3))
    W31 = 4.*x3*x1/(1.-(x3-x1)*(x3-x1))

    G = g_ref + g_ent + W12*G_excess12 + W23*G_excess23 + W31*G_excess31
    expr = G
    
    # Get the derivatives as strings with replacements
    result = partial_derivatives_to_cpp(expr, variables)
    
    # Print the results
    print("Energy:")
    print(f"this->phase_free_energy = {int_to_float(wrap_decimal_numbers(int_to_float(replace_powers_with_intpow(expr))))};\n")
    print("Partial Derivatives:")
    for var, derivative_str in result.items():
        print(f"dfd{var} = {derivative_str};\n")

