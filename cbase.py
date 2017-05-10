#########################################################
#						cBase2							#
#					Syed Abid Hasan						#
#ChemSpider Token: ccab768d-273b-41cb-a63e-2b8ebf2508ad #
#########################################################

#Chemspider queries: http://blog.matt-swain.com/post/16893587098/chemspipy-a-python-wrapper-for-the-chemspider 
#Chemspider queries: http://chemspipy.readthedocs.io/en/latest/guide/gettingstarted.html (new API here)
#rdkit http://www.rdkit.org/docs/api/rdkit.Chem-module.html
#rdkit http://www.rdkit.org/docs/GettingStartedInPython.html
#csv reader http://www.pythonforbeginners.com/systems-programming/using-the-csv-module-in-python/

#Load the CSV Structure file into memory.
# TO--DO: Write a module for converting SDF to CSV, if this is required.
# Should be doable from SQL server though....

#import csv
from cbase_functions import *

#TO--DO : VALIDATE THE FILE (ensure dollar separated, ensure 17 rows or whatever)
#TO--DO : VALIUDATE/LOAD SETTINGS FILE, USERNAME FILE, DICTIONARY

#Login System TO--DO: move this to a function in cbase_functions
valid_name = False		#variable for storing whether it's a valid login name
login_name = ""			#given login name
while valid_name == False:
	if login_name != "":
		#User has supplied an incorrect name before if we are here
		print "Username %s not found." % login_name
	else:
		#first time here
		print "Welcome to cBase v2\nWhat is your name?\n"
	login_name = raw_input("\nPlease enter your username, or type \"e\" to exit > ")
	if login_name.lower() == "e": exit()
	#Open username file
	user_names = open('usernames.txt')		#file stores all usernames separated by "$" character
	user_names_list = user_names.read().split("$")
	for line in user_names_list:
		if login_name == line:
			#we found a matching user to what is in file!
			valid_name = True
			user_names.close()
			break		#exit loop
user_names.close()
	
#todo variable stores what the user wants to do
todo = 0
while todo != constants("exit"):	#TO--DO is this a good way to do things?
	result = []

	if todo != 0:
		#entered something before
		print "\n'%s' is an invalid input! \n\nWhat would you like to do %s?" % (todo, login_name)
	else:
		print "\n\nWelcome %s! What would you like to do? " % login_name
	
	print "1. Search for item by CAS Number"            						###DONE###
	print "2. Search for item by structure"						    		#SEARCHING DONE, CHECK FOR ERROR VALENCE BS
	print "3. Search for item by name or common name"						###DONE###
	print "4. Search for item by chemical formula"                         			#SEARCHING DONE, CHECK FOR ERROR VALENCE BS
	print "5. Search for item by chemical ID"                             						###DONE###
	print "6. Search for item txt, csv file"          	#TO--DO
	print "7. Conduct a common reaction (Mitsonobu, Suzuki, etc.) - not implemented yet"
	print "8. Fix the database (remove quant = 0, unassigned location, structure <=> CAS matching, remove NA cas numbers, fix no mol generation"		#use CAS to SMILES form CACTUS
																					#cactus.nci.nih.gov/chemical/structure and PUBCHEM for NA cas numbers (from name)
	print "9. Preferences (options file edit)"
	print "10. Rebuild dict, fix database"										###DONE###
	print "11. Exit"
	
	todo = ask_input("    >>", "integer", False)

	#OPTION 9
	while todo == constants("housekeeping"):
		rebuild = raw_input("\nAre you sure you want to rebuild the dictionary. This operation may take a while (y/n)?").lower()
		if rebuild == "y":
			dictionary_rebuild()
		else:
			print "Dictionary will not be rebuilt\n\n"				
			todo = 0
		
		fix_cas = raw_input("\nAre you sure you want to fix CAS numbers? All 'NA' CAS numbers will be rebuilt (y/n)").lower()
		if fix_cas == "y":
			cas_rebuild()
		else:
			print "CAS numbers will not be fixed\n\n"
			
		#fix_structures => fix the structure so it doesnt give stupid errors, based on cas number.
	
	#OPTIONS 1, 2, 3, 4, 5
	while todo == constants("search_cas") or todo == constants("search_smiles") or todo == constants("search_name") or todo == constants("search_formula") or todo == constants("search_id"):
		#We are searching by CAS number, Structure/SMILES, or name/common-name
		if todo == constants("search_cas"):
			message = "CAS number"
		elif todo == constants("search_smiles"):
			message = "SMILES string or file"
		elif todo == constants("search_name"):
			message = "Chemical Name"
		elif todo == constants("search_formula"):
			message = "chemical formula"
		elif todo == constants("search_id"):
			message = "Chemical ID"
            
		print "\n\n\nWhich %s would you like to search for?\nOr Type \"e\" to go back to main menu" % message
		input_string = raw_input("%s     >>" % message)
		if input_string.lower() == "e":
			print "\n\n\n----------------------------------------------------"				
			todo = 0		#so we start the main loop fresh (otherwise it'd bitch that you entered wrong input, evne though it was just "1")
		
        #if its CAS, call the CAS searcher, if it's name, get the namesearaher
        
		if todo == constants("search_cas"):
			result = cas_find(input_string)
			#call the CAS Number Searcher, returns a list containing indices of found data
			#returns False if nothing found
		elif todo == constants("search_smiles"):
			result = structure_find(input_string)
		elif todo == constants("search_name"):
			result = name_find(input_string)
		elif todo == constants("search_id"):
			result = id_find(input_string)
		elif todo == constants("search_formula"):
			result = formula_find(input_string)


        #result is either FALSE (nothing found) or a list containing tuples of all found things, along with match score 
		if len(result) > 0:
			#Call function to print relevant results
			find_compound(result, True)		#true means to show the header line, printer determines whether you want results printed to a file
			
			#cpd found, ask what to do with it
			truncated_results = pick_compounds(result)
			
			if truncated_results == "e":
				print "\n\n\n----------------------------------------------------"				
				todo = 0		#so we start the main loop fresh (otherwise it'd bitch that you entered wrong input, evne though it was just "1")
				break
			
			#now let's work with truncated results
			#print "1. Get chemical properties"
			#print "2. GHS Safety data"
			#print "3. Add to shopping cart"
			print "\n" * 5
			print "you picked %s" % truncated_results	
		elif len(result) == 0:
			print "No compounds found with that %s" % message
			
			#TO--DO: this code is here so that when you get to the end of the search, it returns tyo main menu rather than asking again for another search, wbhich is pointless
	#print "\n\n\n----------------------------------------------------"				
	#todo = 0		#so we start the main loop fresh (otherwise it'd bitch that you entered wrong input, evne though it was just "1")		
			
print "Bye, thanks for using cBase!"