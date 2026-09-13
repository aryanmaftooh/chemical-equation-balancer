import re
import numpy as np
from  collections import deque
from scipy.linalg import null_space
from scipy.linalg import qr

def quilifyInputEqu(equation):
    equation_format = r"(?:(\d*(?:[A-Z][a-z]?\d*)+\+)+)?\d*(?:[A-Z][a-z]?\d*)+=(?:(\d*(?:[A-Z][a-z]?\d*)+\+)+)?\d*(?:[A-Z][a-z]?\d*)+"
    prog = re.compile(equation_format)
    return prog.match(equation)

    

def Matrix(num_atom:int, num_molecule:int):
    matrix = []
    for i in range(num_atom):
        matrix.append([])
        for j in range(num_molecule):
            matrix[i].append(0)
    return matrix   
def isatom(atom):
    return atom.isdigit()==False


def findmultiplier(molecule: str)->int:
    multiplier = 1
    findmulti  = re.search(r'[A-Za-z\(\)]?(\d*)', molecule).group(1)
    if findmulti != '':
        multiplier = int(findmulti)
    return multiplier


def findmultiplierParanthesis(molecule: str)->int:
    stack = deque()
    stack.append('(')
    i =  0
    for c in molecule:
        if c == '(':
            stack.append(c)
        elif c == ')':
            stack.pop()
        if len(stack) == 0:
            multiplier = findmultiplier(molecule[i:])
            return multiplier

        i += 1

    try:    
        raise IndexError
    except IndexError:
        raise IndexError(f"Open parantheses (... at {molecule}")


def modifyMatrix(matrix, reactant_molecules, product_molecules, library):
    stack = deque()
    column = 0
    multiplier = 1
    for molecule in reactant_molecules:
        i = 0
        for atom in molecule:
            if atom >= 'a' and atom <= 'z':
                i+=1
                continue
            elif atom == '(':
                t = findmultiplierParanthesis(molecule[i+1:])
                multiplier = multiplier * t
                stack.append(t)
                temp = ''
            elif atom == ')':
                try:
                    multiplier =int(multiplier / stack.pop())
                except IndexError:
                    raise IndexError(f"Open parantheses ...) at {molecule}")
                temp = ''
            elif isatom(atom):
                temp = re.search(fr'{atom}([a-z])?',molecule[i:]).group(0)
                matrix[library[temp]][column]= 1*multiplier
            elif temp != '' and atom.isdigit()==True:
                    matrix[library[temp]][column] = multiplier* int(re.search(fr"{temp}(\d*)", molecule).group(1))
                    temp = ''
            i+=1
        column+=1

    stack1 = deque()
    column =  len(reactant_molecules)-1
    multiplier = 1
    for molecule in product_molecules:
        temp = ''
        i = 0
        column+=1

        for atom in molecule:
             if atom >= 'a' and atom <= 'z':
                i+=1
                continue
             elif atom == '(':
                t = findmultiplierParanthesis(molecule[i+1:])
                multiplier = multiplier * t
        
                stack1.append(t)
                temp = ''
             elif atom == ')':
                try:
                    multiplier =int(multiplier / stack1.pop())
                except IndexError:
                    raise IndexError(f"Open parantheses ...) at {molecule}")
                temp = ''
             elif isatom(atom):
                 temp = re.search(fr'{atom}([a-z])?',molecule[i:]).group(0)
                 matrix[library[temp]][column]=-1*multiplier
             elif   temp != '' and atom.isdigit()==True:
                matrix[library[temp]][column] =  -1 *multiplier* int(re.search(fr"{temp}(\d*)", molecule).group(1))
                temp = ''
             i+=1
def solveNullSpace(matrix):
        null_vec = null_space(np.array(matrix))
        null_vec = np.absolute(null_vec)
        return np.round(null_vec[:, 0] / np.min(null_vec[null_vec > 0])).astype(int)



def openFile():
    file = open("equations", "r")
    for line in file:
        print(line.strip())
    file.close()


def main(chemequ):
        if quilifyInputEqu(chemequ.replace('(', '').replace(')', '').replace(' ', '')) == None:
            print("I am not able to solve this. Please enter a valid equation")
            print("Problem Occured: ", chemegu)
            return 0
        equation = re.search(r'(.*)=(.*)', chemequ)
        reactant = equation.group(1)
        product = equation.group(2)
        reactant_molecules = re.findall(r'[\(\)A-Za-z0-9]+', reactant) #list
        product_molecules = re.findall(r'[\(\)A-Za-z0-9]+', product) #list
        atoms = list(set(re.findall(r'[A-Z][a-z]?', chemequ)))
        matrix = Matrix(len(atoms), len(product_molecules) + len(reactant_molecules))
        dictionary = {}
        counter = 0
        for i in atoms:
            dictionary[i] = counter
            counter +=1
        modifyMatrix(matrix, reactant_molecules, product_molecules, dictionary)
        vec = solveNullSpace(matrix)
        for i in range(len(reactant_molecules)):
            if i != len(reactant_molecules)-1:
                print(f'{vec[i]}{reactant_molecules[i]} +', end= ' ')
        else:
            print(f'{vec[i]}{reactant_molecules[i]}', end = ' ')
        print(" => ", end= ' ')
        for i in range(len(product_molecules)):
            if i != len(product_molecules)-1:
                print(f'{vec[i+len(reactant_molecules)]}{product_molecules[i]} +', end= ' ')
        else:
            print(f'{vec[i+len(reactant_molecules)]}{product_molecules[i]}')
                
   
   



