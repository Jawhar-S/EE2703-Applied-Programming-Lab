'''
----------------------------------------------------------------------------------------------------------------
 EE2703 - Assignment-2
 Done by Jawhar S (EE20B049)
 Last modified on: 09/02/2022
 Description: This program accepts a netlist file in command-line; parses the netlist file to identify different
    components and their attributes in the circuit; constructs the nodal analysis matrix and the source 
    matrix; solves for the nodal voltages and current through voltage sources for both ac and dc circuits. 
----------------------------------------------------------------------------------------------------------------
'''

from sys import argv, exit
import numpy as np

# Identifiers for the start and end of circuit and AC line 
CIRCUIT='.circuit'
END='.end'
AC='.ac'

# Class for storing various parameters of all elements
class allElements():
    
    def __init__(self, line):
        self.line = line
        self.tokens = self.line.split()
        self.name = elementType(self.tokens[0])
        self.fromNode = self.tokens[1]
        self.toNode = self.tokens[2]
        
        # RLC 
        if self.tokens[0][0] in 'RLC'and len(self.tokens) == 4:
            self.type = 'RLC'
            self.value = float(self.tokens[3])

        # DC Sources
        elif self.tokens[0][0] in 'VI' and len(self.tokens) == 5:
            self.type = 'dc'
            self.value = float(self.tokens[4])

        # CCxS
        elif self.tokens[0][0] in 'HF' and len(self.tokens) == 5:
            self.type = 'CCxS'
            self.v_cont = self.tokens[3]
            self.value = float(self.tokens[4])

        # AC Sources
        elif self.tokens[0][0] in 'VI' and len(self.tokens) == 6:
            self.type = 'ac'
            V = float(self.tokens[4])/2
            phase = float(self.tokens[5])*np.pi/180
            real = V*np.cos(phase)
            imaginary = V*np.sin(phase)
            self.value = complex(real, imaginary)

        # VCxS
        else:
            self.type = 'VCxS'
            self.v_node1 = self.tokens[3]
            self.v_node2 = self.tokens[4]
            self.value = float(self.tokens[5])

# Function for getting Element Names
def elementType(token):

    if token[0] == 'R':
        return 'Resistor'
    elif token[0] == 'L':
        return 'Inductor'
    elif token[0] == 'C':
        return 'Capacitor'
    elif token[0] == 'V':
        return 'Independent Voltage Source'
    elif token[0] == 'I':
        return 'Independent Current Source'
    elif token[0] == 'E':
        return 'Voltage Controlled Voltage Source'
    elif token[0] == 'F':
        return 'Current Controlled Current Source'
    elif token[0] == 'G':
        return 'Voltage Controlled Current Source'
    elif token[0] == 'H':
        return 'Current Controlled Voltage Source'
    

# Function to return dictionary of nodes
def nodeDictionary(cktDefn):

    dictionary = {}
    nodes = [allElements(line).fromNode for line in cktDefn]
    nodes.extend([allElements(line).toNode for line in cktDefn])
    i = 1
    nodes = list(set(nodes))
    for node in nodes:
        if node == 'GND':
            dictionary[node] = 0                        # The GND node is assigned value 0
        else:
            dictionary[node] = i
            i += 1
    return dictionary

# Function to get the corresponding key for a value in the dictionary
def key(dictionary, value):

    for key in dictionary.keys():
        if dictionary[key] == value:
            return key

# Function to make a dictionary for each component of a particular type of element
def makeDictionary(cktDefn, element):

    elementDictionary = {}
    elementNames = [allElements(line).tokens[0] for line in cktDefn if allElements(line).tokens[0][0] == element]
    for i,name in enumerate(elementNames):
        elementDictionary[name] = i
    return elementDictionary

def makeList(cktDefn, element):

    elementSet = [allElements(line) for line in cktDefn if allElements(line).tokens[0][0] == element]
    return elementSet

# Function to get number of nodes; voltage sources; VCV sources; CCV sources
def countNodesAndElements(cktDefn):

    voltageIndependent = [i for i in range(len(cktDefn)) if cktDefn[i].split()[0][0] == 'V']
    voltageDependentVoltage = [i for i in range(len(cktDefn)) if cktDefn[i].split()[0][0] == 'E']
    currentDependentVoltage = [i for i in range(len(cktDefn)) if cktDefn[i].split()[0][0] == 'H']
    VS = len(voltageIndependent)
    N = len(nodeDictionary(cktDefn))
    VCVS = len(voltageDependentVoltage)
    CCVS = len(currentDependentVoltage)
    return N, VS, VCVS, CCVS

# Function for finding the position of a given node
def nodePosition(cktDefn, nodeKey, dictionary):

    x = [(i,j) for i in range(len(cktDefn)) for j in range(len(cktDefn[i].split())) if cktDefn[i].split()[j] in dictionary.keys() if dictionary[cktDefn[i].split()[j]] == nodeKey]
    return x

# Function for forming matrices M and b for a given node
def FormMatrices(cktDefn, frequency, nodeKey, dictionary, voltageDictionary, fullVoltageList, inductorDictionary, VCVSDictionary, CCVSDictionary, M, b):

    x = nodePosition(cktDefn,nodeKey,dictionary)
    N,VS,VCVS,CCVS = countNodesAndElements(cktDefn)
    for i in x:
        element = allElements(cktDefn[i[0]])
        elementName = cktDefn[i[0]].split()[0]
        
        # Resistor
        if elementName[0] == 'R':

            if i[1] == 1:
                adjKey = dictionary[element.toNode]
                M[nodeKey,nodeKey] += 1/(element.value)
                M[nodeKey,adjKey] -= 1/(element.value)
            if i[1] == 2 :
                adjKey = dictionary[element.fromNode]
                M[nodeKey,nodeKey] += 1/(element.value)
                M[nodeKey,adjKey] -= 1/(element.value)

        # Capacitor
        if elementName[0] == 'C':
            if i[1] == 1:
                adjKey = dictionary[element.toNode]
                M[nodeKey,nodeKey] += complex(0, 2*np.pi*frequency*(element.value))
                M[nodeKey,adjKey] -= complex(0, 2*np.pi*frequency*(element.value))
            if i[1] == 2 :
                adjKey = dictionary[element.fromNode]
                M[nodeKey,nodeKey] += complex(0, 2*np.pi*frequency*(element.value))
                M[nodeKey,adjKey] -= complex(0, 2*np.pi*frequency*(element.value))

        # Inductor
        if elementName[0] == 'L':
            try:
                if i[1] == 1:
                    adjKey = dictionary[element.toNode]
                    M[nodeKey,nodeKey] -= complex(0,1/(2*np.pi*frequency*element.value))
                    M[nodeKey,adjKey] += complex(0,1/(2*np.pi*frequency*element.value))
                if i[1] == 2 :
                    adjKey = dictionary[element.fromNode]
                    M[nodeKey,nodeKey] -= complex(0,1/(2*np.pi*frequency*element.value))
                    M[nodeKey,adjKey] += complex(0,1/(2*np.pi*frequency*element.value))

            # In case of DC circuits
            except ZeroDivisionError:
                index = inductorDictionary[elementName]
                if i[1] == 1:
                    M[nodeKey,N+VS+VCVS+CCVS+index] += 1 
                    M[N+VS+VCVS+CCVS+index,nodeKey] -= 1
                    b[N+VS+VCVS+CCVS+index] = 0
                if i[1] == 2:
                    M[nodeKey,N+VS+VCVS+CCVS+index] -= 1
                    M[N+VS+VCVS+CCVS+index,nodeKey] += 1
                    b[N+VS+VCVS+CCVS+index] = 0

        # Independent Voltage Source
        if elementName[0] == 'V':
            index = voltageDictionary[elementName]
            if element.value != 0:
                if i[1]== 1:
                    M[nodeKey,N+index] -= 1
                    M[N+index,nodeKey] -= 1
                    b[N+index] = element.value
                if i[1] == 2 :
                    M[nodeKey,N+index] += 1
                    M[N+index,nodeKey] +=1
                    b[N+index] = element.value

        # Independent Current Source
        if elementName[0] == 'I':
            if i[1]== 1:
                b[nodeKey] -= element.value
            if i[1] == 2 :
                b[nodeKey] += element.value

        # Voltage Controlled Voltage Source
        if elementName[0] == 'E':
            index = VCVSDictionary[elementName]
            if i[1]== 1:
                M[nodeKey,N+VS+index] += 1
                M[N+VS+index,nodeKey] += 1
            if i[1] == 2 :
                M[nodeKey,N+VS+index] -= 1
                M[N+VS+index,nodeKey] -=1
            if i[1] == 3:
                M[N+VS+index,nodeKey] += element.value
            if i[1] == 4:
                M[N+VS+index,nodeKey] -= element.value

        # Current Controlled Current Source
        if elementName[0] == 'F':
            index = voltageDictionary[element.v_cont]
            for j in range(len(fullVoltageList)):
                if element.v_cont == fullVoltageList[j].tokens[0]:
                    controlling_n1 = dictionary[fullVoltageList[j].fromNode]
                    controlling_n2 = dictionary[fullVoltageList[j].toNode]
            if i[1]== 1:
                M[nodeKey,N+index] -= element.value
                M[controlling_n1, N+index] -= 1
                M[N+index, controlling_n1] -= 1
            if i[1] == 2:
                M[nodeKey,N+index] += element.value
                M[controlling_n2, N+index] += 1
                M[N+index, controlling_n2] += 1
        
        # Voltage Controlled Current Source
        if elementName[0] == 'G':
            controlling_n1 = dictionary[element.v_node1]
            controlling_n2 = dictionary[element.v_node2]
            if i[1]== 1:
                M[nodeKey,controlling_n1] += element.value
                M[nodeKey,controlling_n2] -= element.value
            if i[1] == 2 :
                M[nodeKey,controlling_n1] -= element.value
                M[nodeKey, controlling_n2] += element.value

        # Current Controlled Voltage Source
        if elementName[0] == 'H':
            index1 = voltageDictionary[element.v_cont]
            index2 = CCVSDictionary[elementName]
            for j in range(len(fullVoltageList)):
                if element.v_cont == fullVoltageList[j].tokens[0]:
                    controlling_n1 = dictionary[fullVoltageList[j].fromNode]
                    controlling_n2 = dictionary[fullVoltageList[j].toNode] 
            if i[1] == 1:
                M[nodeKey,N+VS+VCVS+index2] -= 1
                M[N+VS+VCVS+index2,nodeKey] -= 1
                M[N+VS+VCVS+index2,N+index1] -= element.value
                M[controlling_n1,N+index1] -= 1
                M[N+index1,controlling_n1] -= 1 
            if i[1] == 2:
                M[nodeKey,N+VS+VCVS+index2] += 1
                M[N+VS+VCVS+index2,nodeKey] += 1
                M[controlling_n2,N+index1] += 1
                M[N+index1,controlling_n2] += 1
            

# Function for getting the frequency and circuit definition
def GetFreqCktDefn(netlist_file):

    if (not netlist_file.endswith(".netlist")):
                raise ValueError
    else:
        with open (netlist_file, "r") as f:
            all_lines = []
            for line in f.readlines():
                all_lines.append(line.split('#')[0].split('\n')[0])
            
            start_loc = all_lines.index(CIRCUIT)
            end_loc = all_lines.index(END)
            freq = 0
            ACList = []
            DCList = []
            
            for i in range(len(all_lines)):
                if AC in all_lines[i]:
                    freq = float(all_lines[i].split()[2])
                    ACList.append(freq)
                if 'dc' in all_lines[i]:
                    DCList.append(i)

            ACSet = set(ACList)
            if len(ACSet) >0 and len(DCList) >0:
                exit("Error: This program cannot solve circuits containing both AC and DC Sources")
            if len(ACSet) >1:
                exit("Error: This program cannot solve circuits containing AC Sources with different frequencies")
            ckt_defn = all_lines[start_loc+1:end_loc]
            return freq, ckt_defn

if __name__  == '__main__':
    
    # Validating The Number Of Arguments
    if len(argv) != 2:
        exit("Invalid number of arguments! Please pass the netlist file as the second argument.")
    else:
        try:
            netlist_file = argv[1]
            try:
                frequency, cktDefn = GetFreqCktDefn(netlist_file)
                dictionary = nodeDictionary(cktDefn)
                voltageDictionary = makeDictionary(cktDefn, 'V')
                fullVoltageList = makeList(cktDefn, 'V')
                inductorDictionary = makeDictionary(cktDefn, 'L')
                VCVSDictionary = makeDictionary(cktDefn, 'E')
                CCVSDictionary = makeDictionary(cktDefn, 'H')
                N,VS,VCVS,CCVS = countNodesAndElements(cktDefn)
                
                # Initializing the matrices M and b
                dim = N+VS+VCVS+CCVS
                M = np.zeros((dim, dim), dtype = complex)
                b = np.zeros(dim, dtype = complex)
                
                # Adding current through the inductors as variables for purely DC circuits
                if frequency == 0:
                    
                    M = np.zeros((dim+len(inductorDictionary), dim+len(inductorDictionary)), dtype = complex)
                    b = np.zeros(dim+len(inductorDictionary), dtype = complex)

                # Updating the matrices M And b for each node
                for i in range(len(dictionary)):
                    FormMatrices(cktDefn, frequency, i, dictionary, voltageDictionary, fullVoltageList, inductorDictionary, VCVSDictionary, CCVSDictionary, M, b)

                # Setting the GND voltage to 0
                M[0] = 0
                M[0, 0] = 1
                b[0] = 0
        
                try:
                    # Solving The Equation Mx = b
                    x = np.linalg.solve(M,b) 

                except Exception:

                    # The Incidence Matrix is Singular
                    exit('The incidence matrix cannot be inverted as it is singular. Please provide a valid circuit definition')
                
                # Printing Voltage At Each Node
                for j in range(1,N):
                    print("The voltage at node {} is {:.3e}".format(key(dictionary, j), x[j]))

                # Printing Current Through Each Independent Voltage Source
                for j in range(VS):
                    print('The current through source {} is {:.3e}'.format(key(voltageDictionary, j), x[N+j]))

                # Printing Current Through Each Voltage Controlled Voltage Source
                for j in range(VCVS):
                    print('The current through source {} is {:.3e}'.format(key(VCVSDictionary, j), x[N+VS+j]))

                # Printing Current Through Each Current Controlled Voltage Source
                for j in range(CCVS):
                    print('The current through source {} is {:.3e}'.format(key(CCVSDictionary, j), x[N+VS+VCVS+j]))

                # Printing Current Through Each Inductor (only in pure DC circuits)
                if frequency == 0:

                    for j in range(len(inductorDictionary)):
                        print('The current through inductor {} is {:.3e}'.format(key(inductorDictionary, j), x[N+VS+VCVS+CCVS+j]))
            
            except ValueError:
                print("Invalid netlist file! Please make sure that you have entered a valid netlist file.")
        except FileNotFoundError:
            print("Given file does not exist! Please check if you have entered the name of the netlist file correctly.")





    