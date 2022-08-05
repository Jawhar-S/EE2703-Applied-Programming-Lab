'''
-----------------------------------------------------------------------------------------------------
 EE2703 - Assignment-1
 Done by Jawhar S (EE20B049)
 Description: This program accepts a netlist file as a command-line argument and traverses the circuit
  definiton from the last element to the first; saves the element, corresponding nodes, value
  and prints each line with words in reverse order
-----------------------------------------------------------------------------------------------------
'''

from sys import argv, exit

# identifiers indicating the start and end of the circuit definition in the netlist file
START = '.circuit'
END = '.end'

# function for extracting the tokens from a line
def extract_tokens(line):

    words = line.split()
    element = words[0]
    node_1 = words[1]
    node_2 = words[2]

    # R, L, C, Independent Sources
    if(element[0] in 'RLCVI'):
        value = words[3]
        return [element, node_1, node_2, value]

    # CCVS, CCCS
    elif(element[0] in 'EG'):
        voltage_source = words[3]
        value = words[4]
        return [element, node_1, node_2, voltage_source, value]

    # VCVS, VCCS
    elif(element[0] in 'HF'):
        voltage_source_node_1 = words[3]
        voltage_source_node_2 = words[4]
        value = words[5]
        return [element, node_1, node_2, voltage_source_node_1, voltage_source_node_2, value]

    else:
        return []

# function which returns the list which contains the lines between .circuit and .end of the netlist file
def get_lines_and_words(netlist_file):
    if (not netlist_file.endswith(".netlist")):
                raise ValueError
    else:
        with open (netlist_file, "r") as f:
            all_lines = []
            for line in f.readlines():
                all_lines.append(line.split('#')[0].split('\n')[0])
            
            start_loc = all_lines.index(START)
            end_loc = all_lines.index(END)
            ckt_defn = all_lines[start_loc+1:end_loc]
            lines_and_words = [extract_tokens(line) for line in ckt_defn]
            return lines_and_words

# function for printing the last line to the first with words in each line reversed 
def print_reversed(lines_and_words):
    lines_and_words = reversed(lines_and_words)
    for x in lines_and_words:
        x = reversed(x)
        for y in x:
            print(y, end=' ')
        print('')
    print('')
    return

if __name__ == "__main__":

    # checking number of command line arguments
    if len(argv)!=2 :
        exit("Invalid number of arguments! Please pass the netlist file as the second argument.")
    else:
        try:
            netlist_file = argv[1]
            try:
                lines_and_words = get_lines_and_words(netlist_file)
                print_reversed(lines_and_words)
            except ValueError:
                print("Invalid netlist file! Please make sure that you have entered a valid netlist file.")
        except FileNotFoundError:
            print("Given file does not exist! Please check if you have entered the name of the netlist file correctly.")




