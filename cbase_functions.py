def constants(parameter):
	if parameter == "casrow": return 2
	if parameter == "idrow": return 4			#id row in cas window
	if parameter == "namerow": return 3			#name row in cas window
	if parameter == "smilesrow": return 1			#name row in cas window

	
	if parameter == "exit": return 11
	if parameter == "housekeeping": return 10			#name row in cas window
	if parameter == "search_id": return 5
	if parameter == "search_formula": return 4
	if parameter == "search_name": return 3
	if parameter == "search_smiles": return 2
	if parameter == "search_cas": return 1			#options in main menu
		
	#TO--DO: generate these parameters dynamically
	#TO--DO: make this into a dictionary, or read from settings file? 

def ask_input(text, type, exit):
	'''This function takes a text to display (text), a datatype ("integer"), and exit code (whether 'e' input
	should exit or be ignored, and returns the integer of input'''
	while True:
		entered = raw_input(text)
		#if exit is requested, then if "e' was entered exit and return
		if exit == True:
			if entered == "e": return "e"
		
		if type == "integer":
			try:
				return int(entered)
			except ValueError:
				continue
			
		
def cas_find(cas_number_input):
	#This function takes the cas number supplied, sees if it's valid, then 
	#searches through the  given database file for that number and returns
	#the results as a list
	
	#Check valid CAS
	import re		#regex module
	match_cas = re.search(r'\d{2,7}-\d{2}-\d', cas_number_input, re.I)
	if match_cas:
		#print "Valid CAS number entered (%s).\n" % cas_number_input

		import csv
		csv_file = open('cbase.dollar', 'rb')
		database = csv.reader(csv_file, delimiter='$')

		search_results = []
		#Search through CAS database
		for rownum, row in enumerate(database):
			#print row
			if rownum != 0:
				#non header row, search through
				## TO--DO : row[2] is hard coded for file. Fix this later
				if cas_number_input in row[constants("casrow")]:		#TO--DO: or use "if row[constants("casrow")] == cas_number_input:	...	"
					match_score = round( float(len(cas_number_input)) / float(len(row[constants("casrow")]))*100, 0)
					search_results.append((rownum, match_score))		#this will return a list containing the index numbers
            
		csv_file.close()
		#TO--DO: ensure that csv file is always closed before function returns shit, otherwise it'll never be run!
		#if len(search_results) > 0:
		return search_results
		#else:
			#nothing was found
		#	return []
	else:
		#return False
		print "Entered CAS Number is malformed. CAS numbers should be [######]-##-#!\n\n"

		
def id_find(id_string):
	#This function takes the ID number supplied, then searches through the
	#given database file for that number and returns results as list
	
	if not id_string:
		print "No chemical ID was supplied!\n\n"
		return []
		
	import csv
	csv_file = open('cbase.dollar', 'rb')
	database = csv.reader(csv_file, delimiter='$')

	search_results = []
	#Search through CAS database
	for rownum, row in enumerate(database):
		if id_string in row[constants("idrow")]:
			match_score = round( float(len(id_string)) / float(len(row[constants("idrow")]))*100, 0)
			search_results.append((rownum, match_score))		#this will return a list containing the index numbers

	return search_results
	csv_file.close()


def name_find(name_string):
	'''Function finds a chemical by name. Searches through names from ChemicalBook and in house database. Expects string input, outputs list containing index number of chemical database'''
	from difflib import SequenceMatcher
	import csv
	
	name_string = str(name_string).lower()
	if (not name_string) or len(name_string) < 4:
		#nothing was given, so fail
		print "Invalid or too short name supplied!\n\n"
		return []	
		
	#TO--DO: check for freshness of dictionary file, and prompt rebuild if needed!
	#search through dictionary
	dict_file = open('namesdictionary.dollar')
	dict_db = csv.reader(dict_file, delimiter='$')
	dict_matches = set()		#this will hold the cas numbers that are matched in the dictionary
	
	for counter, item in enumerate(dict_db):
		if counter < 2: continue	#ignore the first two lines	TO--DO: all of these such things should be something like "if left$(line, 1) = "#" then continue to allow multiple commentable rows"
		
		if name_string in str(item[1]).lower():
			dict_matches.add(item[0])

	#now loop through CBASE and look for 1) the name as entered 	and 	2) if the row matches the CAS number from dict_mathces
	csv_file = open('cbase.dollar', 'rb')
	database = csv.reader(csv_file, delimiter='$')
	
	#this set will store the search results (set forces only unique items in here!)
	search_results = []
	
	#search for compounds in in house database
	for rownum, row in enumerate(database):
		if (name_string in row[constants("namerow")]) or (row[constants("casrow")] in dict_matches):
			#we found a hit
			#TO--DO: ADD LOGIC FOR PROPER MATCH PERCENTAGE:::
			#max (  [SequenceMatcher(None, name_string, row[constants("namerow")]), SequenceMatcher(None, name_string, dict_matches))
			search_results.append((rownum, 0))

	csv_file.close()
	dict_file.close()

	return search_results
	
	
def structure_find(structure_input):
	#import RDKit packages:
	from rdkit import Chem
	from rdkit.Chem import Descriptors
	from rdkit.Chem import Lipinski
	from rdkit.Chem import Draw
	from rdkit.Chem.Draw import MolDrawing
	from collections import defaultdict
	#MolDrawing.elemDict=defaultdict(lambda : (0,0,0)) 			#if you want black atom labels

	#parse the input string, if it is none, then return notrhing
	search_pattern = Chem.MolFromSmiles(str(structure_input))	
	if search_pattern is None:
		print "The search string %s is malformed" % structure_input
		return []
	search_heavy_atoms = search_pattern.GetNumHeavyAtoms()
	
	#open the file and load into memory
	import csv
	csv_file = open('cbase.dollar', 'rb')
	database = csv.reader(csv_file, delimiter='$')
	#load the strucutres and IDs
	search_results = []
	
	for rownum, row in enumerate(database):		
		#header row is useless so ignore it:
		if rownum < 1: continue
		#get the SMILES string from database, iterate if string is empty
		molecule  = row[constants("smilesrow")]
		if molecule == "" or molecule is None: continue
		
		#generate molecule object in RDKit from this molecule
		mol = Chem.MolFromSmiles(molecule)
		
		#critical: if generation of molecule object failed, then continue
		if mol is None:
			continue
		
		#match the input pattern given
		match = mol.GetSubstructMatch(search_pattern)
		
		#add the ID row to list if match successful
		if match:
			mol_heavy_atoms = mol.GetNumHeavyAtoms()
			match_percentage = int (search_heavy_atoms / mol_heavy_atoms * 100)
			#TO--DO: ADD LOGIC FOR PROPER MATCH PERCENTAGE. JUST SUM NUMBER OF ATOMS IN QUERY  and DIVIDE BY NUM ATOMS In PDT 
			search_results.append((rownum, match_percentage))

			
#TO--DO: PROPERTIES OF COMPOUNDS
#		if match:
#			count = count + 1
#			print "True",
#			print Descriptors.MolWt(mol),
#			print Lipinski.NumHDonors(mol),
#			print Lipinski.NumHAcceptors(mol),
#			print Lipinski.NumHeteroatoms(mol),
#			print Lipinski.RingCount(mol)
			#Draw.MolToFile(mol,'test.png',size=(300,300))#, highlightAtoms=match)
			
	csv_file.close()	
	return search_results

		
def find_compound(index_list, show_header):
	'''This function prints data about compounds given its index number(s)'''
	#it expects a list, through which it'll loop and show the data
	#Keyword parameters omit certain things from being printed
		#USER, SMILES, CASNUMBER, NAME, ID, MOLWT, AMTORDERED, AMOUNTUNIT, DESCR, DATERECEIVED, LOCATION, DATEORDERED)
		#TO--DO whichever one of above is passed as FALSE, do not print it. if no parameters poassed, then print everything!
		#TO--DO	work on spacing for this!!
	import csv
	csv_file = open('cbase.dollar', 'rb')
	database = csv.reader(csv_file, delimiter='$')		#open the dollar delimited file
		
	#Search for compounds
	search_hits = 0		#number of hits
	to_print = ""		#what will need to be printed
	for i, row in enumerate(database):
		#loop through the database
		if i == 0 and show_header == True:
			to_print += "  ".join(row) + "\n"
		for k in index_list:
			if i == k[0]:
				search_hits += 1
				to_print += str(search_hits) + ".   ".join(row) + "    " + str(k[1]) + "\n"

#		if i in index_list_rows:
#			search_hits += 1
#			to_print += str(search_hits) + ".   ".join(row) + "    " + index_list_score[] + "\n"
	
	print "%d search results were found:\n\n" % search_hits
	if search_hits > 0:
		print to_print
		write_to_file(to_print)
		print "-------------------------------------------------"
	csv_file.close()


def formula_find(input_formula):
	#############################
	# make dictionary of input  #
	#############################
	#if no input_formula or too short then stop
	if input_formula is None or len(input_formula) < 2:
		print "Not a valid formula - too short"
		return []
     
	#get rid of spaces if they exist
	import re
	if re.search(r"[^a-zA-Z0-9]", input_formula) or re.search(r"[^A-Z]", input_formula[0:1]):
		print "Not a valid formula - contains invalid characters"
		return []

	valid_elements = ['Ac', 'Ag', 'Al', 'Am', 'Ar', 'As', 'At', 'Au','B','Ba','Be','Bh','Bi','Bk','Br','C','Ca','Cd',
	'Ce','Cf','Cl','Cm','Co','Cr','Cs','Cu','D','Db','Ds','Dy','Er','Es','Eu','F','Fe','Fm','Fr','Ga','Gd','Ge','H','He',
	'Hf','Hg','Ho','Hs','I','In','Ir','K','Kr','La','Li','Lr','Lu','Md','Mg','Mn','Mo','Mt','N','Na','Nb','Nd','Ne',
	'Ni','No','Np','O','Os','P','Pa','Pb','Pd','Pm','Po','Pr','Pt','Pu','Ra','Rb','Re','Rf','Rg','Rh','Rn','Ru','S',
	'Sb','Sc','Se','Sg','Si','Sm','Sn','Sr','Ta','Tb','Tc','Te','Th','Ti','Tl','Tm','U','Uu','Uu','Uu','Uu','Uu','Uu',
	'Uu','V','W','Xe','Y','Yb','Zn','Zr']
    
	#find all elements using regex:      returns a list with tuples:   [('S' : 1), ('H' : 2), ('X' : 5)]
	input_reg = re.findall(r"([A-Z][a-z]?)(\d{0,2})", input_formula)
	#create a dictionary for storing results
	dict_input = dict()
	#loop through list, and find the elements from tuples
	for item in input_reg:
		print item
	#item will be ('O', '4') tuple
	#number of the current element (eg, O4 ===> '4'). If nothing, then '1'
		if item[1] != '':
			num = int(item[1])
		else:
			num = 1
        
        #check if this element is in list of valid elements
		if (item[0] in valid_elements) == False:
			print "Element %s is not valid" % item[0]
			return []

	#add to ductionary if its not there, otherwise update the count
		if item[0] in dict_input:
			dict_input[item[0]] += num
		else:
			dict_input[item[0]] = num
    ##########################################
    # loop through db, make dict of entries  #
    ##########################################

	#open and loop through file
	import csv
	from rdkit import Chem
	from rdkit.Chem.rdMolDescriptors import CalcMolFormula

	csv_file = open('cbase.dollar', 'rb')
	database = csv.reader(csv_file, delimiter='$')

    #to store the resutls taht will be passed back
	search_results = []
	for rownum, row in enumerate(database):
		if rownum < 1: continue
        
		molecule  = str(row[constants("smilesrow")])
		if molecule == "" or molecule is None: continue
        
		curr_mol = Chem.MolFromSmiles(molecule)
		if curr_mol is None: continue

		#get formula
		curr_formula = str(CalcMolFormula(curr_mol))
		
		#make dictionary of formula
		curr_reg = re.findall(r"([A-Z][a-z]?)(\d{0,2})", curr_formula)
		#create a dictionary for storing results        
		dict_curr = dict()
        
		for curr_item in curr_reg:
			num = curr_item[1]
			if num =='': num = 1

			if curr_item[0] in dict_curr:
				dict_curr[curr_item[0]] += int(num)
			else:
				dict_curr[curr_item[0]] = int(num)
        
		#MATCH THE CUYRR DICTIOANRY WITH THE INPUT DICTIONARY
		#Attempt to match the two. If lengths are different then obviously there is no match
		if len(dict_curr) != len(dict_input): continue

		#loop through the dictionary. key loops through the keys of the dictionary
		match = True
		
		for key in dict_input: 
			if not(str(key) in dict_curr):
			#print "no key in curr dict"
				match = False
				break
		if match == False: continue
		
		if dict_input[key] != dict_curr[key]:
			match = False
			break
		if match == False: continue
		#TO--DO: ADD LOGIC FOR PROPER MATCH PERCENTAGE
		if match == True: search_results.append((rownum, 1))

	csv_file.close()
	return search_results
	

def dictionary_rebuild():
	'''This function rebuild the dictionary REAL'''
	import csv		#for reading CSV file
	import re		#for finding valid CAS numbers
	import sys
	
	csv.field_size_limit(2147483647)		#TO--DO: this is a windows specific issue!
	
	cbase_file = open('cbase.dollar')
	cbase_db = csv.reader(cbase_file, delimiter='$')
	
	dict_file = open('namesdictionarymaster.dollar')
	dict_db = csv.reader(dict_file, delimiter='$')
	
	#read cbase into array
	cbase_list = []
	for counter, item in enumerate(cbase_db):
		if counter <1: continue			#if first row then ignore - its the title row
		curr_cas = str(item[constants("casrow")])
		match_cas = re.search(r'\d{2,7}-\d{2}-\d', curr_cas, re.I)
		if match_cas: cbase_list.append(curr_cas)
	
	#loop through major (read only) dictionary
	matches = []
	version = ""
	for counter, item in enumerate(dict_db):
		if counter < 2: 
			version += str(item[0]) + "\n"		#version number, etc lines (first two are useless lines)
			continue
			
		if len(item) != 2:
			print "current row is messed up %s : %s" % (counter, item)
			raw_input(" > ")
			continue
		curr_cas = str(item[0])
		curr_synonym = str(item[1])
		
		sys.stdout.write("\r currently on line %s, reading %s" % (counter, curr_cas))
		sys.stdout.flush()
		
		#check if in cbase DB
		if curr_cas in cbase_list:
			matches.append( curr_cas + "$" + curr_synonym)
	

	synonyms_file = open('namesdictionary.dollar', 'w')
	#write version number
	synonyms_file.write(version[:-2])		#-2 to prevent the extra new line
	for item in matches:
		synonyms_file.write(str(item) + "\n")
	

	synonyms_file.close()
	cbase_file.close()
	dict_file.close()

def cas_rebuild():
	'''this function aims to look at cbase.dollar and remove all CAS = 'na' and replace them with proper cas numbers'''
	#Load modules that are needed
	import csv		#for reading CSV file
	import re		#for finding valid CAS numbers
	from difflib import SequenceMatcher		#Match text percentages
	import sys
	import wikipedia
	
	csv.field_size_limit(2147483647)		#TO--DO: this needs to be sorted, but it gives an error when large CSVs are loaded if it's omitted

	master_file = open('namesdictionarymaster.dollar')
	master_db = csv.reader(master_file, delimiter='$')
	mast_dict = {}
	
	for counter, x in enumerate(master_db):
		if counter < 2: continue		#first rows are comments, TO--DO: rahter than hard code, do      if  row[0:1] == "#" then continue
		mast_dict[	str(x[1])	] = x[0]

	#this file stores known fixes, so they dont have to be re-read from the master database or manually entered
	knownfix_file = open('knowndictfixes.dollar')
	knownfixes = knownfix_file.readlines()
	known_fixes_dict = {}
	for item in knownfixes:
		if item[0:1] == "#": continue
		line = item.split('$')
		if len(line) == 2:
			known_fixes_dict[line[0]] = line[1].rstrip()
	knownfix_file.close()
	
	#just for fun keep count of what kind of fixes we have
	fixed_count = [0, 0, 0, 0]		#each is a different kind of fix
	
	#loop over cbase	
	cbase_file = open('cbase.dollar')
	cbase_db = csv.reader(cbase_file, delimiter='$')
	recent_fixes = {}
	
	
	for counter, row in enumerate(cbase_db):
		if counter<1: continue
		sys.stdout.write("\rCurrently checking compound %s" % (counter))
		sys.stdout.flush()
		
		curr_cas = str(row[constants("casrow")])
		curr_name = str(row[constants("namerow")])
		if not(curr_cas) and not(curr_name):
			continue
		
	
		match_cas = re.search(r'\d{2,7}-\d{2}-\d', curr_cas, re.I)
		if not(match_cas):			
			knownfix_file = open('knowndictfixes.dollar', 'a')
			
			#see if in dictionary, otherwise loop through
			print "\n\nCurrently on row %s with invalid CAS number '%s' and name '%s'" % (counter, curr_cas, curr_name)
			if curr_name in known_fixes_dict:
				print "cBase name '%s' and CAS '%s' matched in KNOWN FIXES FILE to CAS '%s'" % (curr_name, curr_cas, known_fixes_dict[curr_name])
				fixed_count[0] += 1
			elif curr_name in recent_fixes:
				print "cBase name '%s' and CAS '%s' matched in RECENT FIXES to Master CAS '%s'" % (curr_name, curr_cas, recent_fixes[curr_name])
				fixed_count[1] += 1
			elif curr_name in mast_dict:
				recent_fixes[curr_name] = mast_dict[curr_name]
				#new_known_fixes.append(str(curr_name) + "$" + mast_dict[curr_name])
				print "cBase name '%s' and CAS '%s' matched in MASTER DICTIONARY to CAS '%s'" % (curr_name, curr_cas, mast_dict[curr_name])
				fixed_count[2] += 1
			else:
				#exact match was not found. Do three things
				# 1 - use sequence matcher to loop through master dictionary
				# 2 - ask Wikipedia for CAS number
				# 3 - ask user for CAS number
				best_cas = ""
				#Searching Wikipedia
				wiki_search = wikipedia.search(curr_name)
				if len(wiki_search) > 0:
					print "Found these matches:   " + str(wiki_search)
					wiki_pg = wikipedia.WikipediaPage(wikipedia.search(curr_name)[0])
					try:
						wiki_cas = re.search(r'\d{2,7}-\d{2}-\d', wiki_pg.html(), re.I)
					except:
						pass
					#wiki_syn = re.search(r'other names.*</div>', wiki_pg.html(), re.I)
					if wiki_cas:
						try:
							print "Page title: " + str(wiki_pg.title) + "     CAS:" + str(wiki_cas.group(0))
							best_cas = wiki_cas.group(0)
						except:
							pass
					else:
						print "No CAS found ... " + str(wiki_pg.url)
				
				if len(best_cas) == 0:		#cas wasnt found in Wikipedia
					print "Wikipedia Searching was not fruitful..."
					#Search Sequence Matcher
					print "Attempting to thoroughly search the master dictionary..."
					count = 0
					largest = [0, "", ""]			#closest match found, it is Match Percentage, CAS number, name that matched
					for chem_names in mast_dict:
						count += 1
						if abs(len(chem_names) - len(curr_name)) > 5: continue			#obviously this is going to be no match, no need to do computationally complex stuff
						if largest[1] == mast_dict[chem_names]: continue				#this CAS number already exists, so no need to redo it
					
						ratio = SequenceMatcher(None, chem_names, curr_name).ratio()
						if ratio > largest[0]:
							largest = [ratio, mast_dict[chem_names], chem_names] 
						sys.stdout.write("\rLooping through master, on %s / %s, sequence similarity is %s" % (count, len(mast_dict), ratio))
						sys.stdout.flush()
					print "\nClosest match found (with match percentage of %s) is %s with CAS number %s" % (largest[0], largest[2], largest[1])
					best_cas = largest[1]
			
				while True:
					todo = raw_input("[y] to accept above CAS '%s', or type your own > " % best_cas).lower()
					#todo = "y"
					if todo == "y":
						recent_fixes[curr_name] = str(best_cas) 		#largest[1]
						#new_known_fixes.append(str(curr_name) + "$" + largest[1])
						fixed_count[2] += 1
						break
					elif re.search(r'\d{2,7}-\d{2}-\d', todo, re.I):
						#accept the entered CAS
						recent_fixes[curr_name] = todo
						#new_known_fixes.append(str(curr_name) + "$" + todo)
						print "%s will be used for name %s" % (todo, curr_name)
						fixed_count[3] += 1
						break
				
				line = str(curr_name) + '$' + str(recent_fixes[curr_name])
				knownfix_file.write(line + "\n")
				knownfix_file.close()
				
	#write the known fixes file for next time
		
	print "Fixes from file %d" % fixed_count[0]
	print "Repeated fixes %d" % fixed_count[1] 
	print "Fixes read from master file %d" % fixed_count[2]
	print "Custom entered fixes %d" % fixed_count[3] 
	
	#TO--DO: (1) write the actual cbase file
	
	knownfix_file.close()
	cbase_file.close()
	master_file.close()


def write_to_file(print_text):
	import time
#TO--DO: have an option in here to read a preferecen for whether to write a log file or not....

	#This function writes the given file to text
	#TO--DO: relative path here. Perhaps use settings file or constants function above
	file_write_path = "C:/Python27/Learn/CBASE/Logs/results - %s.txt" % time.strftime("%d-%m-%Y %H-%M-%S")
	file_write = open(file_write_path, 'w')
	file_write.write(print_text)
	file_write.close()
	print "File was written successfully to", file_write_path


def pick_compounds(result):
	'''This function asks for further profiling on already found compounds.'''
	'''Recives a list produced by cas_find, id_find, etc.'''
	#for regex matching of input
	import re
	
	#the message to be displayed, depending on how many things are in the result set
	if len(result) == 1:
		compound_number_msg = "(1)"
	else:
		compound_number_msg = "(1 - %d)" % len(result)

	#Truncated results will hold whatever is selected
	truncated_results = []
	
	while True:
		#keeps track of whether we have a valid range
		valid_range = True
		
		range_input = raw_input("Enter a compound or range of compounds (separated by comma/hyphens), or 'e' for exit, or 'a' for all compounds: %s" % compound_number_msg)
		
		#if we enter "a" for every compound or "e" for exit the return, or nothing entered
		if range_input.lower() == "e": 
			return "e"
		elif range_input.lower() == "a":
			return result
		elif range_input == "":
			print "Nothing was entered!"
			continue
			
		#check for string errors
		#		  non hyphen/comma/numerics	    consecutive non-digit characters	starting with non-digit				ending with non-digits
		if re.search(r"[^0-9,-]", range_input) or re.search(r"[^0-9]{2,}", range_input) or re.search(r"^[^0-9]", range_input) or re.search(r"[^0-9]$", range_input):
			#Invalid entry. Restart the loop for another input. Here, valid_range is only used for completion sake, it isnt really used for anything
			valid_range = False
			print "Invalid range '%s' was found." % range_input
			continue

		#split by commas to get groups of digits
		range_groups = range_input.split(",")
		#now loop through these groups
		for groups in range_groups:
			#either a group is a digit (e.g. '2', '66', '22', etc.) or it's a range '2-6', 6-10', etc.
			if "-" in groups:
				#further split the group into the two ranges
				group_split = groups.split("-")
				
				#check for validity of range
				#			is the lower bound within range		upper bound within range						#lower is less than upper one (or theyre equal)
				if not(int(group_split[0]) in xrange(1, len(result) + 1)) or not(int(group_split[1]) in xrange(1, len(result) + 1)) or int(group_split[0]) >= int(group_split[1]):
					#invalid so break out of this loop (the valid_range = false will force while loop to recycle)
					print "Invalid range entered '%s'" % groups
					valid_range = False
					break
				else:
					#now the range is valid, so lets loop through it
					for i in range(int(group_split[0]), int(group_split[1]) + 1):
						#if this result doesn't already exist, then add it, otherwise complain that it's there
						if not(result[i - 1] in truncated_results):
							truncated_results.append(result[i - 1])
						else:
							print "Selection %s in range '%s' is already chosen!" % (i, groups)
			else:
				#just a single compound, check if it's within the range of result set
				if not (int(groups) in xrange(1, len(result) + 1)):
					print "Selection '%s' is out of range" % groups
					valid_range = False
				else:
					#now the range is valid, lets add it if it doenst already exist
					if not(result[int(groups) - 1] in truncated_results):
						truncated_results.append(result[int(groups) - 1])
					else:
						print "Selection '%s' is already chosen!" % int(groups)
		
		#if valid_range is true, then everyhthig upo there worked out, so lets exit the loop
		if valid_range == True: break
	
	#if result set exists, then return it, otherwise nothing chosen
	if len(truncated_results) > 0:
		return truncated_results
	else:
		print "Nothing valid chosen."
		return False	
	
	
	
	
	
	
	
	
	
	
	
	
