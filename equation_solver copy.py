import tkinter as tk
from tkinter import ttk, messagebox
from tkinter import scrolledtext
from tkinter import *
from sympy import *
import sympy as sym
from sympy.parsing.sympy_parser import (parse_expr, standard_transformations, implicit_multiplication_application)
from sympy import sin, cos, tan, exp, log


# Create a function to solve equation using the bisection method


def bisection(xL, xU, error, max_iterations, equation, rounding_digits):
    """Calculate Xr and Error values."""
    iteration_counter = 0
    
    if max_iterations == 0:
        return False
    
    solution = []
    index = 1
    
    xR = (xL + xU) / 2
    
    # 1# Iteration
    iteration_counter += 1
    if (f(equation, xL) * f(equation, xR)) < 0:
        xU = xR
    else:
        xL = xR
    
    solution.append((index, round(xL, rounding_digits), round(f(equation, xL), rounding_digits), round(xU, rounding_digits), round(f(equation, xU), rounding_digits), round(xR, rounding_digits), round(f(equation, xR), rounding_digits), "---"))
    index += 1
    
    if max_iterations <= iteration_counter:
        return solution
            
    while(True):
        iteration_counter += 1
        xR_Old = xR
        xR = (xL + xU) / 2
        
        solution.append((index, round(xL, rounding_digits), round(f(equation, xL), rounding_digits), round(xU, rounding_digits), round(f(equation, xU), rounding_digits), round(xR, rounding_digits), round(f(equation, xR), rounding_digits), round((abs((xR - xR_Old) / xR) * 100), rounding_digits)) )
        index+=1
        
        if ((abs((xR - xR_Old) / xR) * 100) <= error ):
            return solution;

        if (f(equation, xL) * f(equation, xR)) < 0:
            xU = xR;
        else:
            xL = xR;
        
        if max_iterations <= iteration_counter:
            return solution
        
def false_position(xL, xU, error, max_iterations, equation, rounding_digits):
    """Calculate Xr and Error values."""
    
    iteration_counter = 0
    
    if max_iterations == 0:
        return False
    
    solution = []
    index = 1
    
    xR = xU - ((f(equation, xU) * (xL - xU)) / (f(equation, xL) - f(equation, xU)))
    
    # 1st Iteration
    iteration_counter += 1
    
    if f(equation, xL) * f(equation, xR) < 0.0:
        xU = xR
    else:
        xL = xR
        
    solution.append((
        index,
        round(xL, rounding_digits),
        round(f(equation, xL), rounding_digits),
        round(xU, rounding_digits),
        round(f(equation, xU), rounding_digits),
        round(xR, rounding_digits),
        round(f(equation, xR), rounding_digits),
        "---"
    ))
    index += 1    
    
    if max_iterations <= iteration_counter:
        return solution
            
    while True:
        iteration_counter += 1
        xR_old = xR
        xR = xU - ((f(equation, xU) * (xL - xU)) / (f(equation, xL) - f(equation, xU)))
        eps = abs((xR - xR_old) / xR) * 100
        
        solution.append((
            index,
            round(xL, rounding_digits),
            round(f(equation, xL), rounding_digits),
            round(xU, rounding_digits),
            round(f(equation, xU), rounding_digits),
            round(xR, rounding_digits),
            round(f(equation, xR), rounding_digits),
            str(round(eps, 3)) + "%"
        ))
        index += 1
        
        if abs((xR - xR_old) / xR) * 100 <= error:
            return solution

        if f(equation, xL) * f(equation, xR) < 0.0:
            xU = xR
        else:
            xL = xR
            
        if max_iterations <= iteration_counter:
            return solution
        
        
        

def simplefpoint(x0, error, equation, max_iterations, rounding_digits):
    x = Symbol('x')
    g = equation_converter(equation)
    g_prime = g.diff(x)
    if not g or not g_prime:
        return False
    g_lambda = lambdify(x, g.evalf())
    g_prime_lambda = lambdify(x, g_prime.evalf())
    i = 0
    xi = x0 # Set xi to x0 for the first iteration
    rows = []
    while i < max_iterations:
        i += 1
        x_new = g_lambda(xi)
        error_i = abs((x_new - xi) / x_new)
        row = [i, round(x_new, rounding_digits), round(g_lambda(x_new), rounding_digits), round(error_i, rounding_digits)]
        rows.append(row)
        if error_i < error:
            break
        xi = x_new # Update xi to x_new for the next iteration
    return rows
def newton(xo, error, equation, max_iterations, rounding_digits):
    x = sym.symbols('x')
    f = equation
    f_diff = sym.diff(f, x)
    
    iteration_counter = 0

    if max_iterations == 0:
        return False

    solution = []
    index = 1

    # Iteration #1
    iteration_counter += 1
    xi = xo
    solution.append((index, round(xi, rounding_digits), round(f.evalf(subs={x: xi}), rounding_digits), round(f_diff.evalf(subs={x: xi}), rounding_digits), "---"))
    index += 1

    if max_iterations <= iteration_counter:
        return True

    # Iteration #2
    iteration_counter += 1
    xiplus1 = xi - (f.evalf(subs={x: xi}) / f_diff.evalf(subs={x: xi}))
    eps = abs((xiplus1 - xi) / xiplus1) * 100
    solution.append((index, round(xiplus1, rounding_digits), round(f.evalf(subs={x: xiplus1}), rounding_digits), round(f_diff.evalf(subs={x: xiplus1}), rounding_digits), str(round(eps, rounding_digits)) + "%"))
    index += 1

    if eps <= error:
        return solution

    if max_iterations <= iteration_counter:
        return True

    while eps > error:
        iteration_counter += 1
        xi = xiplus1
        xiplus1 = xi - (f.evalf(subs={x: xi}) / f_diff.evalf(subs={x: xi}))
        eps = abs((xiplus1 - xi) / xiplus1) * 100
        solution.append((index, round(xiplus1, rounding_digits), round(f.evalf(subs={x: xiplus1}), rounding_digits), round(f_diff.evalf(subs={x: xiplus1}), rounding_digits), str(round(eps, rounding_digits)) + "%"))
        index += 1

        if eps <= error:
            return solution

        if max_iterations <= iteration_counter:
            return True
        
        
        
        



# Define the secant method function

from sympy import Symbol  
# Define the secant method function
def secant(xi_minus1, xi, error, max_iterations, rounding_digits, equation):
    iteration_counter = 0
    if max_iterations == 0:
        return False

    solution = []
    index = 1

    while True:
        iteration_counter += 1
        eps = abs((xi - xi_minus1) / xi) * 100

        if iteration_counter == 1:
            solution.append((index, round(xi_minus1, rounding_digits), f(equation, xi_minus1), round(xi, rounding_digits), f(equation, xi), str(round(eps, rounding_digits))+"%"))

        xi_old = xi
        xi = xi - ((f(equation, xi) * (xi_minus1 - xi)) / (f(equation, xi_minus1) - f(equation, xi)))
        xi_minus1 = xi_old

        index += 1
        solution.append((index, round(xi_minus1, rounding_digits), f(equation, xi_minus1), round(xi, rounding_digits), f(equation, xi), str(round(eps, rounding_digits))+"%"))

        if eps <= error:
            return solution

        if iteration_counter >= max_iterations:
            return solution
            
def f(equation, value):
    x = Symbol('x')
    equation = equation.subs('sin', sin) # Substitute sin(x) with sin(x)
    equation = equation.subs('cos', cos)
    equation = equation.subs('tan', tan)
    equation = equation.subs('exp', exp)
    equation = equation.subs('log', log) # Substitute log(x) with log(x)
    parsed_eq = equation_converter(equation)
    return round(parsed_eq.evalf(subs={x: value}), 3)


# Create a function to convert the equation to a format that can be evaluated by SymPy
def equation_converter(eq):
    transformations = (standard_transformations + (implicit_multiplication_application,))
    eq = eq.replace('^', '**')  # Replace ^ with **
    try:
        parsed_eq = parse_expr(eq, transformations=transformations)
        parsed_eq = parsed_eq.subs({'sin': sin, 'cos': cos, 'tan': tan, 'exp': exp, 'log': log}) # Substitute sin(x) with sin(x), cos(x) with cos(x), etc.
        return parsed_eq
    except:
        return False
    
    

# Create a function to check if the function has a solution or not for Bisection and False Position methods
def checkFunctionHasSolutionOrNot(xL, xU, equation):
    if (f(equation, xL) * f(equation, xU) < 0):
        return True
    else:
        return False

def run_program():
    equation = equation_entry.get()
    method = method_var.get()
    error = float(error_entry.get())
    max_iterations = int(max_iterations_entry.get())
    rounding_digits = int(rounding_digits_entry.get())  # prompt user for rounding digits

    # Check if the equation is valid
    parsed_equation = equation_converter(equation)
    if not parsed_equation:
        messagebox.showerror("Error", "Invalid equation")
        return
    
    if method == "simple_fixed_point":
        x_value = float(x_value_entry.get())
        g_equation_str = str(parsed_equation) + ' + x'
        solution = simplefpoint(x_value, error, g_equation_str, max_iterations, rounding_digits)  # pass rounding_digits to function call
        simplefp_table.delete(*simplefp_table.get_children())
        for row in solution:
            simplefp_table.insert("", "end", values=row)
        simplefp_frame.pack()
        table_frame.pack_forget()
    elif method == "newton":
        xo = float(x_value_entry.get()) # Use float instead of str
        solution = newton(xo, error, parsed_equation, max_iterations, rounding_digits)  # pass rounding_digits to function call
        newton_table.delete(*newton_table.get_children())
        for row in solution:
            newton_table.insert("", "end", values=row)
        table_frame.pack()
        simplefp_frame.pack_forget()
    else:
        xL = float(xL_entry.get())
        xU = float(xU_entry.get())

        # Check if the function has a solution within the given interval
        if not checkFunctionHasSolutionOrNot(xL, xU, parsed_equation):
            messagebox.showerror("Error", "Function has no solution within the given interval")
            return

        # Run the selected method
        if method == "bisection":
            solution = bisection(xL, xU, error, max_iterations, parsed_equation, rounding_digits)  # pass rounding_digits to function call
            bisection_table.delete(*bisection_table.get_children())
            for row in solution:
                bisection_table.insert("", "end", values=row)
            table_frame.pack()
            simplefp_frame.pack_forget()
        elif method == "false_position":
            solution = false_position(xL, xU, error, max_iterations, parsed_equation, rounding_digits)  # pass rounding_digits to function call
            false_position_table.delete(*false_position_table.get_children())
            for row in solution:
                false_position_table.insert("", "end", values=row)
            table_frame.pack()
            simplefp_frame.pack_forget()
        elif method == "secant":
            xi_minus1 = float(xi_minus1_entry.get())
            xi = float(xi_entry.get())
            solution = secant(xi_minus1, xi, error, max_iterations, rounding_digits) # pass rounding_digits to function call
            secant_table.delete(*secant_table.get_children())
            for row in solution:
                secant_table.insert("", "end", values=row)
            table_frame.pack()
            simplefp_frame.pack_forget() 
root = tk.Tk()
root.geometry("800x600")
root.title("Equation Solver")

# Create a frame for the equation and method widgets
equation_frame = tk.Frame(root)
equation_frame.pack(fill=tk.X, padx=10, pady=10)

# Create the equation label and entry widgets
equation_label = tk.Label(equation_frame, text="Equation:")
equation_label.pack(side=tk.LEFT, padx=5)
equation_entry = tk.Entry(equation_frame, width=50)
equation_entry.pack(side=tk.LEFT, padx=5)

# Create the method label and dropdown widgets
method_label = tk.Label(equation_frame, text="Method:")
method_label.pack(side=tk.LEFT, padx=5)
method_var = tk.StringVar()
method_dropdown = tk.OptionMenu(equation_frame, method_var, "bisection", "false_position","simple_fixed_point","newton","secant")
method_dropdown.pack(side=tk.LEFT, padx=5)

# Create a frame for the parameter widgets
parameter_frame = tk.Frame(root)
parameter_frame.pack(fill=tk.X, padx=10, pady=10)

# Create the xL label and entry widgets
xL_label = tk.Label(parameter_frame, text="xL:")
xL_label.pack(side=tk.LEFT, padx=5)
xL_entry = tk.Entry(parameter_frame, width=10)
xL_entry.pack(side=tk.LEFT, padx=5)

# Create the xU label and entry widgets
xU_label = tk.Label(parameter_frame, text="xU:")
xU_label.pack(side=tk.LEFT, padx=5)
xU_entry = tk.Entry(parameter_frame, width=10)
xU_entry.pack(side=tk.LEFT, padx=5)

# Create the x0 label and entry widgets
x_value_label = tk.Label(parameter_frame, text="x0:")
x_value_label.pack(side=tk.LEFT, padx=5)
x_value_entry = tk.Entry(parameter_frame, width=10)
x_value_entry.pack(side=tk.LEFT, padx=5)

# Create the xi-1 label and entry widgets
xi_minus1_label = tk.Label(parameter_frame, text="xi-1:")
xi_minus1_label.pack(side=tk.LEFT, padx=5)
xi_minus1_entry = tk.Entry(parameter_frame, width=10)
xi_minus1_entry.pack(side=tk.LEFT, padx=5)

# Create the xi label and entry widgets
xi_label = tk.Label(parameter_frame, text="xi:")
xi_label.pack(side=tk.LEFT, padx=5)
xi_entry = tk.Entry(parameter_frame, width=10)
xi_entry.pack(side=tk.LEFT, padx=5)

# Create the error label and entry widgets
error_label = tk.Label(parameter_frame, text="Error:")
error_label.pack(side=tk.LEFT, padx=5)
error_entry = tk.Entry(parameter_frame, width=10)
error_entry.pack(side=tk.LEFT, padx=5)

# Create the max iterations label and entry widgets
max_iterations_label = tk.Label(parameter_frame, text="Max Iterations:")
max_iterations_label.pack(side=tk.LEFT, padx=5)
max_iterations_entry = tk.Entry(parameter_frame, width=10)
max_iterations_entry.pack(side=tk.LEFT, padx=5)

# Create the rounding digits label and entry widgets (for the simple fixed point method)
simplefp_frame = tk.Frame(root)
simplefp_frame.pack(fill=tk.X, padx=10, pady=10)
simplefp_label = tk.Label(simplefp_frame, text="Simple Fixed Point Method")
simplefp_label.pack(side=tk.LEFT, padx=5)
rounding_digits_label = tk.Label(simplefp_frame, text="Rounding Digits:")
rounding_digits_label.pack(side=tk.LEFT, padx=5)
rounding_digits_entry = tk.Entry(simplefp_frame, width=10)
rounding_digits_entry.pack(side=tk.LEFT, padx=5)

# Create the run button
run_button = tk.Button(root, text="Run", command=run_program)
run_button.pack(side=tk.BOTTOM, padx=10, pady=10)

# Create a frame for the table
table_frame = tk.Frame(root)
table_frame.pack(fill=tk.BOTH, padx=10, pady=10)

# Create the bisection table
bisection_table = ttk.Treeview(table_frame, yscrollcommand=1)
bisection_table['columns'] = ('i', 'xL', 'f(xL)', 'xU', 'f(xU)', 'xR', 'f(xR)', 'eps')
bisection_table.column("#0", width=0, stretch=tk.NO)
bisection_table.column("i", anchor=tk.CENTER, width=80)
bisection_table.column("xL", anchor=tk.CENTER, width=80)
bisection_table.column("f(xL)", anchor=tk.CENTER, width=80)
bisection_table.column("xU", anchor=tk.CENTER, width=80)
bisection_table.column("f(xU)", anchor=tk.CENTER, width=80)
bisection_table.column("xR", anchor=tk.CENTER, width=80)
bisection_table.column("f(xR)", anchor=tk.CENTER, width=80)
bisection_table.column("eps", anchor=tk.CENTER, width=80)
bisection_table.heading("#0", text="", anchor=tk.CENTER)
bisection_table.heading("i", text="Id", anchor=tk.CENTER)
bisection_table.heading("xL", text="xL", anchor=tk.CENTER)
bisection_table.heading("f(xL)", text="f(xL)", anchor=tk.CENTER)
bisection_table.heading("xU", text="xU", anchor=tk.CENTER)
bisection_table.heading("f(xU)", text="f(xU)", anchor=tk.CENTER)
bisection_table.heading("xR", text="xR", anchor=tk.CENTER)
bisection_table.heading("f(xR)", text="f(xR)", anchor=tk.CENTER)
bisection_table.heading("eps", text="eps", anchor=tk.CENTER)
scrollbar = Scrollbar(table_frame, orient=tk.VERTICAL, command=bisection_table.yview)
bisection_table.config(yscrollcommand=scrollbar.set)
scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
bisection_table.pack(fill=tk.BOTH, expand=1)


false_position_table = ttk.Treeview(table_frame, yscrollcommand=1)
false_position_table['columns']= ('i', 'xL', 'f(xL)', 'xU', 'f(xU)', 'xR', 'f(xR)', 'eps')
false_position_table.column("#0", width=0, stretch=NO)
false_position_table.column("i", anchor=CENTER, width=80)
false_position_table.column("xL", anchor=CENTER, width=80)
false_position_table.column("f(xL)", anchor=CENTER, width=80)
false_position_table.column("xU", anchor=CENTER, width=80)
false_position_table.column("f(xU)", anchor=CENTER, width=80)
false_position_table.column("xR", anchor=CENTER, width=80)
false_position_table.column("f(xR)", anchor=CENTER, width=80)
false_position_table.column("eps", anchor=CENTER, width=80)
false_position_table.heading("#0", text="", anchor=CENTER)
false_position_table.heading("i", text="Id", anchor=CENTER)
false_position_table.heading("xL", text="xL", anchor=CENTER)
false_position_table.heading("f(xL)", text="f(xL)", anchor=CENTER)
false_position_table.heading("xU", text="xU", anchor=CENTER)
false_position_table.heading("f(xU)", text="f(xU)", anchor=CENTER)
false_position_table.heading("xR", text="xR", anchor=CENTER)
false_position_table.heading("f(xR)", text="f(xR)", anchor=CENTER)
false_position_table.heading("eps", text="eps", anchor=CENTER)
false_position_table.pack_forget()
scrollbar = Scrollbar(table_frame, orient=tk.VERTICAL, command=bisection_table.yview)
false_position_table.config(yscrollcommand=scrollbar.set)
scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
false_position_table.pack(fill=tk.BOTH, expand=1)


simplefp_table = ttk.Treeview(simplefp_frame, yscrollcommand=1)
simplefp_table['columns'] = ('i', 'xi', 'f(xi)', 'eps')
simplefp_table.column("#0", width=0, stretch=NO)
simplefp_table.column("i", anchor=CENTER, width=80)
simplefp_table.column("xi", anchor=CENTER, width=80)
simplefp_table.column("f(xi)", anchor=CENTER, width=80)
simplefp_table.column("eps", anchor=CENTER, width=80)
simplefp_table.heading("#0", text="", anchor=CENTER)
simplefp_table.heading("i", text="Id", anchor=CENTER)
simplefp_table.heading("xi", text="xi", anchor=CENTER)
simplefp_table.heading("f(xi)", text="f(xi)", anchor=CENTER)
simplefp_table.heading("eps", text="eps", anchor=CENTER)
scrollbar = Scrollbar(root, orient=VERTICAL, command=simplefp_table.yview)
simplefp_table.config(yscrollcommand=scrollbar.set)
scrollbar.pack(side=RIGHT, fill=Y)
simplefp_table.pack()

secant_table = ttk.Treeview(table_frame, yscrollcommand=1)
secant_table['columns'] = ('i', 'xi_minus1', 'f(xi_minus1)', "xi",'f(xi)',  'eps') 
secant_table.column("#0", width=0, stretch=NO)
secant_table.column("i", anchor=CENTER, width=80)
secant_table.column("xi_minus1", anchor=CENTER, width=80)
secant_table.column("f(xi_minus1)", anchor=CENTER, width=80)
secant_table.column("xi", anchor=CENTER, width=80)
secant_table.column("f(xi)", anchor=CENTER, width=80)
secant_table.column("eps", anchor=CENTER, width=80)
secant_table.heading("#0", text="", anchor=CENTER)
secant_table.heading("i", text="Id", anchor=CENTER)
secant_table.heading("xi_minus1", text="xi", anchor=CENTER)
secant_table.heading("f(xi_minus1)", text="f(xi)", anchor=CENTER)
secant_table.heading("xi", text="f'(xi)", anchor=CENTER)
secant_table.column("f(xi)", anchor=CENTER, width=80)
secant_table.heading("eps", text="eps", anchor=CENTER)
scrollbar = Scrollbar(root, orient=VERTICAL, command=bisection_table.yview)
secant_table.config(yscrollcommand=scrollbar.set)
scrollbar.pack(side=RIGHT, fill=Y)
secant_table.pack()
table_frame.pack()

newton_table = ttk.Treeview(table_frame, yscrollcommand=1)
newton_table['columns'] =('i', 'xi', 'f(xi)', "f'(xi)",  'eps')
newton_table.column("#0", width=0, stretch=NO)
newton_table.column("i", anchor=CENTER, width=80)
newton_table.column("xi", anchor=CENTER, width=80)
newton_table.column("f(xi)", anchor=CENTER, width=80)
newton_table.column("f'(xi)", anchor=CENTER, width=80)
newton_table.column("eps", anchor=CENTER, width=80)
newton_table.heading("#0", text="", anchor=CENTER)
newton_table.heading("i", text="Id", anchor=CENTER)
newton_table.heading("xi", text="xi", anchor=CENTER)
newton_table.heading("f(xi)", text="f(xi)", anchor=CENTER)
newton_table.heading("f'(xi)", text="f'(xi)", anchor=CENTER)
newton_table.heading("eps", text="eps", anchor=CENTER)
scrollbar = Scrollbar(root, orient=VERTICAL, command=bisection_table.yview)
newton_table.config(yscrollcommand=scrollbar.set)
scrollbar.pack(side=RIGHT, fill=Y)
newton_table.pack()
table_frame.pack()

def method_changed(*args):
    """Hide/show the table and input widgets based on the selected method."""
    if method_var.get() == "bisection":
        bisection_table.pack()
        false_position_table.pack_forget()
        simplefp_table.pack_forget()
        newton_table.pack_forget()
        xL_entry.pack()
        xU_entry.pack()
        x_value_entry.pack_forget()
    elif method_var.get() == "false_position":
        bisection_table.pack_forget()
        false_position_table.pack()
        simplefp_table.pack_forget()
        newton_table.pack_forget()
        secant_table.pack_forget()
        xL_entry.pack()
        xU_entry.pack()
        xi_minus1_entry.pack_forget() 
        xi_entry.pack_forget()
        x_value_entry.pack_forget()
    elif method_var.get() == "newton":
        bisection_table.pack_forget()
        false_position_table.pack_forget()
        simplefp_table.pack_forget()
        newton_table.pack()
        xL_entry.pack_forget()
        xU_entry.pack_forget()
        x_value_entry.pack()
        xi_minus1_entry.pack_forget() 
        xi_entry.pack_forget()
    elif method_var.get() == "simple_fixed_point":
        bisection_table.pack_forget()
        false_position_table.pack_forget()
        simplefp_table.pack()
        xL_entry.pack_forget()
        xU_entry.pack_forget()
        x_value_entry.pack()
    elif method_var.get() == "secant":
        bisection_table.pack_forget()
        false_position_table.pack_forget()
        newton_table.pack_forget()
        simplefp_table.pack_forget()
        secant_table.pack()
        xL_entry.pack_forget()
        xU_entry.pack_forget()
        x_value_entry.pack_forget()
        xi_minus1_entry.pack()
        xi_entry.pack()  
method_var.trace('w', method_changed)

root.mainloop()

