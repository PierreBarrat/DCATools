"""
SUPPORT LIBRARY
Generic-purpose helper functions.

Sections:
  1. PDB parsing
  2. Error logging
  3. Downloading
  4. Format editing
"""

### 0. Required packages

import os
import sys
import time
import codecs
import datetime
import pickle
import sqlite3
import multiprocessing
import numpy as np


### 1. PDB parsing

# Amino acid properties
# WARNING: it needs the hardcoded file aminoacid_stats.csv!
def read_aa_stats(aa_stats_filename='aminoacid_stats.csv'):
	"""
	Returns the conversion dictionaries from 3 to 1-letter aminoacid notation
	and viceversa, and a dictionary with the number of heavy atoms per residue.
	Warning: it reads a hardcoded file which must be in its same volder
	"""
	# The same file also contains pKa, buried area, etc
	this_dir = sys.path[0]
	aa_stats_file = open(aa_stats_filename, 'r')
	text = aa_stats_file.read().split('\n')
	aa_stats_file.close()

	three2one = {}
	one2three = {}
	nHH = {}
	for line in text:
		if not line:
			continue
		if line.startswith('#'):
			continue
		fields = line.split(',')
		oneletter = fields[0].strip()
		threeletters = fields[1].strip()
		nC = int(fields[11].strip())
		nN = int(fields[13].strip())
		nO = int(fields[14].strip())-1	# Also the final OH is counted
		nS = int(fields[15].strip())

		three2one[threeletters] = oneletter
		one2three[oneletter] = threeletters
		nHH[oneletter] = nC + nN + nO + nS

	return (three2one, one2three, nHH)


# PDB parsing
def pdb_parser_stub(locations, pdb_filename, pdbname='', pkl='', thr_log_status='ERROR'):
	"""
	PDB parser.
	Warning: this is just a stub! This parser is not good for a general PDB file
	It also makes use of the "locations" dictionary, which is available only if you use
	the environment of these python scripts collections
	It records each atom with coordinates (its chain, residue type, whether the residue is
	completely described or some atoms are missing - except for the hydrogens and deuteriums)
	If there is no chain name, the assigned chain name is "_"
	"""
	pdb_file = open(pdb_filename, 'r')
	text = pdb_file.read().split('\n')
	pdb_file.close()

	three2one, one2three, nHH = read_aa_stats()

	pdb_info = {'AUTHOR_CHAINS' : set(), 'DESCRIBED_RESIDUES' : {}}
	iHH = 0
	for line in text:
		if not line:
			continue
		if line.startswith('HEADER') and not pdbname:
			pdbname = line.split()[-1].lower()
		if 'ENDMDL' in line:
			break
		if not line.startswith('ATOM'):
			continue
		altloc = line[16]
		resname = line[17:20]
		chain = line[21]
		if chain == ' ':
			chain = '_'
		resid = line[22:26].strip()
		atname = line[11:16].strip()
		coords = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
		if altloc != ' ' or resname not in three2one or line[26] != ' ' or atname == 'OXT':
			continue
		pdb_info['AUTHOR_CHAINS'].add(line[21])
		if (chain, resid) not in pdb_info['DESCRIBED_RESIDUES']:
			pdb_info['DESCRIBED_RESIDUES'][(chain, resid)] = {'RES_TYPE' : resname, 'HAS_FULL_COORDS' : False, 'COORDS' : {}}
			pdb_info['DESCRIBED_RESIDUES'][(chain, resid)]['COORDS'] = {}
			iHH = 0
		if atname[0] != 'H' and atname[0] != 'D':	# D is deuterium in hydrogen diffraction experiments
			iHH += 1
		if iHH == nHH[three2one[resname]]:
#			print_log(('NOTICE', pdb_parser_stub.__name__, '{0}_{1} res {2} has full coordinates'.format(pdbname, chain, resid)), thr_log_status)
			pdb_info['DESCRIBED_RESIDUES'][(chain, resid)]['HAS_FULL_COORDS'] = True
		elif iHH > nHH[three2one[resname]]:
			print_log(('ERROR', pdb_parser_stub.__name__, 'While analyzing file {0}\nToo many heavy atoms in chain {1} res {2} {3}\n{4} > {5}'.format(pdb_filename, chain, resname, resid, iHH, nHH[three2one[resname]])), thr_log_status)
		pdb_info['DESCRIBED_RESIDUES'][(chain, resid)]['COORDS'][atname] = coords
#		print_log(('NOTICE', pdb_parser_stub.__name__, '{0}_{1} res {2} atoms with coords {3} expected {4}'.format(pdbname, chain, resid, iHH, nHH[three2one[resname]])), thr_log_status)

	if pkl:
		print_log(('NOTICE', pdb_parser_stub.__name__, 'Writing analysis of PDB {0} (PDB file {1}) on {2}'.format(pdbname, pdb_filename, pkl)), thr_log_status)
		pickle.dump(pdb_info, open(locations['tmp']['pdbs'] + pdbname + '.pkl', 'wb'))

	return pdb_info


def process_complex_header_parts(process_text, pdb_header, key_set=False, thr_log_status='ERROR'):
	"""
	General PDB parser subroutine that processes some standard header parts
	that can vary a lot between different PDB files.
	It records them in a dictionary with the following keys:
	MOLECULES: description of each molecule appearing in the PDB (its ID, the
	name of the chain if any, the organism from which it was taken, its TAXID,
	its strain)
	RESOLUTION: the resolution in angstrom, in case of crystallographic structures
	BIOLOGICAL_ASSEMBLIES: all listed biological assemblies listed in order of appearence,
	with author/software indication, chains where it is applied, and biomatrix transformations
	(coming soon) MISSING_RESIDUES: a list of residues not resolved in the structure
	(coming soon) MISSING_ATOMS: a list of atoms not described in the PDB
	(coming soon) SEQUENCE: the complete sequence of the protein, as it was expressed by the
	organism
	"""
	process_keys = set(['COMPND', 'SOURCE', 'REMARK   2', 'REMARK   3', 'REMARK 350', 'REMARK 465', 'REMARK 470', 'SEQRES'])
	if key_set:
		return process_keys

	# TO DO!
	pdb_header['MOLECULES'] = []
	molecule = {}
	text = process_text['COMPND'].split('\n')
	for line in text:
		if not line:
			continue
		if 'MOL_ID:' in line:
			pdb_header['MOLECULES'].append({})
			molecule = pdb_header['MOLECULES'][-1]
			molecule['ID'] = int(line[line.index('MOL_ID:') + 7:-1].strip())
		elif 'MOLECULE:' in line:
			molecule['NAME'] = line[line.index('MOLECULE:') + 9:-1].strip()
		elif 'CHAIN:' in line:
			molecule['AUTHOR_CHAINS'] = set([x.strip() for x in line[line.index('CHAIN:') + 6:-1].strip().split(',')])

	text = process_text['SOURCE'].split('\n')
	for line in text:
		if not line:
			continue
		if 'MOL_ID:' in line:
			mol_id = int(line[line.index('MOL_ID:') + 7:-1].strip())
			n_entry = -1
			for nmol, mol in enumerate(pdb_header['MOLECULES']):
				if mol['ID'] == mol_id:
					n_entry = nmol
					break
			if n_entry == -1:
				pdb_header['MOLECULES'].append({})
				molecule = pdb_header['MOLECULES'][-1]
				molecule['ID'] = mol_id
			else:
				molecule = pdb_header['MOLECULES'][n_entry]
		elif 'ORGANISM_SCIENTIFIC:' in line:
			molecule['ORGANISM'] = line[line.index('ORGANISM_SCIENTIFIC:') + 20:-1].strip()
		elif 'ORGANISM_TAXID:' in line:
			molecule['TAXID'] = line[line.index('ORGANISM_TAXID:') + 15:-1].strip()
		elif 'STRAIN:' in line:
			molecule['STRAIN'] = line[line.index('STRAIN:') + 7:-1].strip()
		
	text = process_text['REMARK   2'].split('\n')
	for line in text:
		if not line:
			continue
		if 'RESOLUTION.' in line and 'ANGSTROMS.' in line:
			pdb_header['RESOLUTION'] = float(line.split()[-1])

	pdb_header['BIOLOGICAL_ASSEMBLIES'] = []
	text = process_text['REMARK   350'].split('\n')
	for line in text:
		if not line:
			continue
		if 'BIOMOLECULE:' in line:
			pdb_header['BIOLOGICAL_ASSEMBLIES'].append({'AUTHOR' : 'NONE', 'SOFTWARE' : 'NONE', 'APPLY_TO' : set(), 'BIOMATRICES' : []})
			apply_lines = False
		elif 'AUTHOR DETERMINED BIOLOGICAL UNIT:' in line:
			pdb_header['BIOLOGICAL_ASSEMBLIES'][-1]['AUTHOR'] = line[line.index('BIOLOGICAL UNIT:') + 16:].strip()
		elif 'SOFTWARE DETERMINED QUATERNARY STRUCTURE:' in line:
			pdb_header['BIOLOGICAL_ASSEMBLIES'][-1]['SOFTWARE'] = line[line.index('QUATERNARY STRUCTURE:') + 21:].strip()
		elif 'APPLY THE FOLLOWING TO CHAINS:' in line:
			pdb_header['BIOLOGICAL_ASSEMBLIES'][-1]['APPLY_TO'] = set([x.strip() for x in line[line.index('QUATERNARY STRUCTURE:') + 21:].strip().split(',') if x.strip()])
			apply_lines = True
		elif 'BIOMT' in line and apply_lines:
			fields = line.split()
			if 'BIOMT1' in line:
				biomatrix = {'ROT' : np.zeros((3,3)), 'TRANSL' : np.zeros((3))}
				biomatrix['ROT'][0][0] = float(fields[4])
				biomatrix['ROT'][0][1] = float(fields[5])
				biomatrix['ROT'][0][2] = float(fields[6])
			elif 'BIOMT2' in line:
				biomatrix['ROT'][1][0] = float(fields[4])
				biomatrix['ROT'][1][1] = float(fields[5])
				biomatrix['ROT'][1][2] = float(fields[6])
			elif 'BIOMT3' in line:
				biomatrix['ROT'][2][0] = float(fields[4])
				biomatrix['ROT'][2][1] = float(fields[5])
				biomatrix['ROT'][2][2] = float(fields[6])
				pdb_header['BIOLOGICAL_ASSEMBLIES'][-1]['BIOMATRICES'].append(biomatrix)

	# Missing residues (to retrieve complete sequence when SEQRES is not there)
	text = process_text['REMARK   465'].split('\n')
	for line in text:
		if not line:
			continue

	# Missing atoms (to check which residues are not complete)
	text = process_text['REMARK   470'].split('\n')
	for line in text:
		if not line:
			continue

	# To retrieve the complete sequence	
	text = process_text['REMARK   465'].split('\n')
	for line in text:
		if not line:
			continue

	return pdb_header


def PDB_parser(pdb_filename, thr_log_status='ERROR'):
	"""
	General PDB parser, complete with all header information and
	keywords from the PDB files, compiled in an orderly way.
	It does not need anything more than the PDB file.

	"""
	pdb_info = {'HEADER' : {}, 'COORDS' : {}}
	header_keys = {'HEADER_TITLE' : ('HEADER', (10,50)), 
	               'DATE' : ('HEADER', (50, 58)),
	               'PDB_ID' : ('HEADER', (62, 66)),
	               'PUB_TITLE' : ('TITLE', (10,80)),
	               'EXP_TECHNIQUE' : ('EXPDTA', (10, 80)),
	               'AUTHOR' : ('AUTHOR', (10, 80))}
	header_keywords = set()
	for k in header_keys:
		header_keywords.add(header_keys[k][0])
	process_keys = process_complex_header_parts('', {}, key_set=True)
	coord_6c_keys = set(['ATOM', 'TER', 'HETATM', 'MASTER', 'MODEL', 'CONECT', 'ANISOU', 'ENDMDL', 'END'])

	atom_line = {'LINETYPE' : 'ATOM', 'COORDS' : np.zeros(3), 'RESNAME' : 'XXX', 'RESID' : 0, 'ATID' : 0,  'ATNAME' : 'XXX', 'CHAIN' : '0', 'ALTLOC' : '0'}
	hetatom_line = {'LINETYPE' : 'HETATOM', 'COORDS' : np.zeros(3), 'RESNAME' : 'XXX', 'RESID' : 0, 'ATID' : 0,  'ATNAME' : 'XXX', 'CHAIN' : '0', 'ALTLOC' : '0'}
	anisou_line = {'LINETYPE' : 'ANISOU', 'REF_ATID' : 0, 'U' : np.zeros((3,3))}
	ter_line = {'LINETYPE' : 'TER'}
	model_line = {'LINETYPE' : 'MODEL', 'MODEL_ID' : 'X'}
	endmdl_line = {'LINETYPE' : 'ENDMDL'}
	end_line = {'LINETYPE' : 'END'}
	master_line = {'LINETYPE' : 'MASTER'}
	

	# 1. Parses file
	pdb_file = open(pdb_filename, 'r')
	text = pdb_file.read().split("\n")
	pdb_file.close()

	process_text = {}
	coords_per_line = {}
	is_header_section = True
	pdb_header = pdb_info['HEADER']
	body_lines = {}
	models_present = False
	for nl, line in enumerate(text):
		if not line:
			continue
		if is_header_section:
			for coord_6c_key in coord_6c_keys:
				if line.startswith(coord_6c_key):
					# Before going on, process all the header keys that were recorded
					pdb_header = process_complex_header_parts(process_text, pdb_header)
					is_header_section = False
			if is_header_section:
				# Parse simple information (info on retrieval contained in header_keys)
				header_keyword = line[:6].strip()
				if header_keyword in header_keywords:
					for k in header_keys:
						if header_keys[k][0] == header_keyword:
							if k not in pdb_header:
								pdb_header[k] = ''
							else:
								pdb_header[k] += ' '
							pdb_header[k] += line[header_keys[k][1][0]:header_keys[k][1][1]].strip()
				for process_key in process_keys:
					if line.startswith(process_key):
						if process_key not in process_text:
							process_text[process_key] = ''
						process_text[process_key] += line + '\n'
				continue
		if line.startswith('ATOM') or line.startswith('HETATM') or line.startswith('ANISOU'):
			if line.startswith('ATOM'):
				body_lines[nl] = atom_line
			elif line.startswith('HETATM'):
				body_lines[nl] = hetatm_line
			else:
				body_lines[nl] = anisou_line
			body_lines[nl]['ATOM_NAME'] = line[11:16].strip()
			body_lines[nl]['ALTLOC'] = line[16]
			body_lines[nl]['RESIDUE_NAME'] = line[17:20]
			body_lines[nl]['CHAIN_NAME'] = line[21]
			body_lines[nl]['RESIDUE_ID'] = int(line[22:26].strip())
			if line.startswith('ANISOU'):
				u00 = float(line[28:35])/10000
				u11 = float(line[36:42])/10000
				u22 = float(line[42:49])/10000
				u01 = float(line[49:56])/10000
				u02 = float(line[56:63])/10000
				u12 = float(line[63:70])/10000
				body_lines[nl]['U'] = np.array([[u00, u01, u02], [u01, u11, u12], [u02, u12, u22]])
			else:
				body_lines[nl]['COORDS'] = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
				body_lines[nl]['OCCUPANCY'] = line[54:60]	# Before transforming in float, check if it really is occupancy
				body_lines[nl]['B_FACTOR'] = line[60:66]	# Before transforming in float, check if it really is B factor
			body_lines[nl]['SEGMENT_ID'] = line[72:76].strip()
			body_lines[nl]['ELEMENT'] = line[76:68].strip()
			body_lines[nl]['CHARGE'] = float(line[78:80])
			if body_lines[nl]['CHAIN_NAME'] == ' ':
				body_lines[nl]['CHAIN_NAME'] = '_'
		elif line.startswith('TER'):
			body_lines[nl] = ter_line		
		elif line.startswith('MASTER'):
			body_lines[nl] = master_line
		elif line.startswith('END'):
			body_lines[nl] = end_line
		elif line.startswith('MODEL'):
			body_lines[nl] = model_line
			fields = line.split()
			if len(fields) == 2:
				body_lines[nl]['MODEL_ID'] = line.split()[1]
				models_present = True
			else:
				print("MODEL LINE WITHOUT A MODEL NUMBER")
		elif line.startswith('CONECT'):
			body_lines[nl] = conect_line
			fields = []
			for pos in range(11, 61, 5):
				fields.append(int(line[pos: pos+5].strip()))
			for fi in range(1, len(fields)):
				body_lines[nl]['CONNECTION'] = (int(fields[0]), int(fields[fi]))
			
		elif line.startswith('ENDMDL'):
			body_lines[nl] = endmdl_line

	# PROCESS COORDINATES
	model_id = ['X', '0']
	model_lines = {}
	for nl in list(sorted(body_lines.keys())):
		if body_lines[nl]['LINETYPE'] == 'MODEL':
			model_id[0] = body_lines[nl]['MODEL_ID']
			model_id[1] = '0'
			if model_id not in model_lines:
				model_lines[model_id] = []
		elif body_lines[nl]['LINETYPE'] == 'ATOM':
			if body_lines[nl]['ALTLOC'] != ' ':
				model_id[1] = body_lines[nl]['ALTLOC']
				if model_id not in model_lines:
					model_lines[model_id] = []
			else:
				model_id[1] = '0'
			model_lines[model_id].append(nl)
		
		
		
		# FIRST LEVEL OF pdb_info['COORDS'] MUST BE MODEL. MAKE 1 MODEL FOR EACH ALTLOC CODE TOO.
		# MODEL NAME: #MODEL_ALTLOC. NO MODEL: X. NO ALTLOC: 0. EX: pdb_info['COORDS']['X_0'], pdb_info['COORDS']['2_A']

	res_ids = {}
	b_factor = {}
	b_norm = {}
	seqres = {}
	applied_chains = []
	biomolecule_occurrences = 0
	tech_list = ['NMR', 'X-RAY', 'THEORETICAL', 'ELECTRON']
	print_log(this_name, "Struct {0}".format(struct))




	for line in text:

		# Technique 'TECHNIQUE'
		if line[0:6] == 'EXPDTA':
			for tech in tech_list:
				if tech in line:
					PDB_dict['TECHNIQUE'] = tech
			if 'TECHNIQUE' not in PDB_dict:
				PDB_dict['TECHNIQUE'] = 'OTHER'
			print_log(this_name, "Struct {0}    Technique {1}".format(struct, PDB_dict['TECHNIQUE']))
		# Resolution 'RESOLUTION'
		elif line[0:10] == 'REMARK   2' and 'RESOLUTION' in line:
			fields = line.split()
			for nf in range(len(fields)):
				if 'ANGSTROM' in fields[nf]:
					if fields[nf-1] != "NULL":
						PDB_dict['RESOLUTION'] = float(fields[nf-1])
					elif not header_from:
						print_log(this_name, "WARNING: Null resolution for structure {0}".format(new_path))
						PDB_dict['RESOLUTION'] = 9999
		# SEQRES
		elif line[0:6] == 'SEQRES':
			fields = line.split()
			ch_name = fields[2]
			if ch_name not in seqres:
				seqres[ch_name] = ''
			seqres[ch_name] += ''.join([from3to1(x) for x in fields[4:]])
		# Free R value 'RFACTOR'
		elif line[0:10] == 'REMARK   3' and 'FREE R VALUE' in line and 'ERROR' not in line and 'SET' not in line and (line.split()[3] == 'FREE' or line.split()[3] == 'BIN'):
			try:
				PDB_dict['RFACTOR'] = float(line.split()[-1])
			except ValueError:
				PDB_dict['RFACTOR'] = 'NULL'
		# (RedIDs, ResNames) - BFactors and BFactorNorms
		elif line[0:4] == 'ATOM':
			if not line[21]:
				raise NameError("ERROR: There is an ATOM without chain name: {0}".format(line))
			ch_name = line[21]
			if ch_name not in res_ids:
				res_ids[ch_name] = []
				b_factor[ch_name] = 0
				b_norm[ch_name] = 0
			if ch_name in res_ids and (not res_ids[ch_name] or  res_ids[ch_name][-1][0] != int(line[22:26])):
				res_ids[ch_name].append((int(line[22:26]), line[17:20]))
			if line[60:66]:
				b_factor[ch_name] += float(line[60:66])
				b_norm[ch_name] += 1
		# Title 'TITLE'
		elif line[0:5] == 'TITLE':
			if 'TITLE' not in PDB_dict:
				PDB_dict['TITLE'] = line[10:].rstrip()
			else:
				PDB_dict['TITLE'] += line[10:].rstrip()
		# Biomolecule 'BIOLOGICAL_UNIT' and biomatrix 'BIOMATRIX'
		elif 'REMARK 350' in line:
			print_log(this_name, "{0}   {1}".format(struct, line))
			if 'BIOMOLECULE:' in line:
				biomolecule_occurrences += 1
			if biomolecule_occurrences > 1:
				continue	
			if 'APPLY THE FOLLOWING TO CHAINS:' in line or 'AND CHAINS:' in line:
				if 'APPLY THE FOLLOWING TO CHAINS:' in line:
					applied_chains = []
				if 'BIOLOGICAL_UNIT' not in PDB_dict:
					PDB_dict['BIOLOGICAL_UNIT'] = set()
				fields = line.split()
				are_chains = False
				for field in fields:
					if are_chains:
						PDB_dict['BIOLOGICAL_UNIT'].add(field[0])
						applied_chains.append(field[0])
					if field == 'CHAINS:':
						are_chains = True
				print("APPLIED TO CHAINS", applied_chains)
			elif 'BIOMT1' in line or 'BIOMT2' in line or 'BIOMT3' in line:
				if not 'BIOMATRIX' in PDB_dict:
					PDB_dict['BIOMATRIX'] = {}
				fields = line.split()
				if not int(fields[3])-1 in PDB_dict['BIOMATRIX']:
					PDB_dict['BIOMATRIX'][int(fields[3])-1] = {}
				Mline = int(fields[2][-1])-1
				PDB_dict['BIOMATRIX'][int(fields[3])-1][Mline] = [float(fields[4]), float(fields[5]), float(fields[6]), float(fields[7])]
				if not applied_chains:
					continue
				PDB_dict['BIOMATRIX'][int(fields[3])-1]['APPLIED_TO'] = applied_chains 
		elif line[0:6] == 'ENDMDL':
			break	# Only the first model is read


	# 2. Fills the keys that have not been found
	#   Null technique 'TECHNIQUE'
	if 'TECHNIQUE' in PDB_dict:
		if PDB_dict['TECHNIQUE'] == 'NMR' and not 'BIOLOGICAL_UNIT' in PDB_dict:
			PDB_dict['BIOLOGICAL_UNIT'] = set()
			for chain in sorted(list(res_ids.keys())):
				PDB_dict['BIOLOGICAL_UNIT'].add(chain)
	else:
		PDB_dict['TECHNIQUE'] = 'UNKNOWN'
	#   Null resolution 'RESOLUTION'
	if 'RESOLUTION' not in PDB_dict and PDB_dict['TECHNIQUE'] != 'THEORETICAL' and PDB_dict['TECHNIQUE'] != 'NMR' and not header_from:
		print_log(this_name, "WARNING: Null resolution for structure {0}".format(new_path))
		PDB_dict['RESOLUTION'] = 9999
	#   Null free R value 'RFACTOR'
	if 'RFACTOR' not in PDB_dict:
		PDB_dict['RFACTOR'] = 'NULL'
	#   Null title 'TITLE':
	if 'TITLE' not in PDB_dict:
		PDB_dict['TITLE'] = ''


	# 3. Fills chainwise data: 'CHAINS', 'AVG_BFACTOR', 'NRES', 'RESIDS', 'RESNAMES'
	PDB_dict['CHAINS'] = set()
	for chain in sorted(list(res_ids.keys())):
		PDB_dict['CHAINS'].add(chain)
		if b_factor[chain] == 0:
			b_factor[chain] = 'NULL'
		else:
			b_factor[chain] = b_factor[chain] / b_norm[chain]
		if chain in seqres:
			sqrs = seqres[chain]
		else:
			sqrs = ''
		PDB_dict[chain] = {'AVG_BFACTOR' : b_factor[chain],
		                   'NRES' : len(res_ids[chain]),
		                   'RESIDS' : [x[0] for x in res_ids[chain]],
		                   'RESNAMES' : [from3to1(x[1]) for x in res_ids[chain]],
		                   'SEQRES' : sqrs}
		print('SEQRES', chain, sqrs)


	# 4. Checks on biological unit 'BIOLOGICAL_UNIT'
	#   If it contains chains not present in 'CHAINS', it eliminates them from 'BIOLOGICAL_UNIT'
	if 'BIOLOGICAL_UNIT' in PDB_dict:
		check_biounit = set()
		for chain in PDB_dict['BIOLOGICAL_UNIT']:
			if chain in PDB_dict['CHAINS']:
				check_biounit.add(chain)
		PDB_dict['BIOLOGICAL_UNIT'] = check_biounit


	if not PDB_dict['CHAINS']:
		print("Something is very wrong with the format of file {0}".format(new_path))
		if pdbfrom:
			print("Since the structure is from an external database, the structure will be scheduled for deletion in OPM_data")
			if hkeys:
				return {}, header_keys
			else:
				return {}
		else:
			print("Since the structure was processed internally, this is an EncoMPASS error and must be reported")
			raise NameError("Format error")

	#   If the biological unit is not present, or if it is present but the biomatrix is not, then it resets the biological unit to 'CHAINS'
	#   and the biomatrix to the identity (applied on all chains)
	if 'BIOLOGICAL_UNIT' not in PDB_dict or ('BIOLOGICAL_UNIT' in PDB_dict and not PDB_dict['BIOLOGICAL_UNIT']) or 'BIOMATRIX' not in PDB_dict:
		if 'BIOLOGICAL_UNIT' in PDB_dict and 'BIOMATRIX' not in PDB_dict:
			print("BIOLOGICAL UNIT PRESENT BUT NOT BIOMATRIX", PDB_dict['BIOLOGICAL_UNIT'])
		PDB_dict['BIOMATRIX'] = {}
		PDB_dict['BIOMATRIX'][0] = {}
		PDB_dict['BIOMATRIX'][0][0] = [1., 0., 0., 0.]
		PDB_dict['BIOMATRIX'][0][1] = [0., 1., 0., 0.]
		PDB_dict['BIOMATRIX'][0][2] = [0., 0., 1., 0.]
		PDB_dict['BIOMATRIX'][0]['APPLIED_TO'] = set(PDB_dict['CHAINS'])
		PDB_dict['BIOLOGICAL_UNIT'] = set(PDB_dict['CHAINS'])
	if pdbfrom and pdbfrom != 'PDB':
		PDB_dict['BIOLOGICAL_UNIT'] = set(PDB_dict['CHAINS'])

	# Eliminates biiomatrices referring to eliminated chains, and corrects APPLIED_TOs
	togo = []
	for n in PDB_dict['BIOMATRIX']:
		new_applied = set()
		for c in PDB_dict['BIOMATRIX'][n]['APPLIED_TO']:
			if c in PDB_dict['CHAINS']:
				new_applied.add(c)
		if new_applied:
			PDB_dict['BIOMATRIX'][n]['APPLIED_TO'] = new_applied
		else:
			print('All chains to which this matrix was applied were eliminated',  PDB_dict['BIOMATRIX'][n])
			togo.append(n)
	for n in togo:
		del PDB_dict['BIOMATRIX'][n]
			
	#    If the biomatrix is present but its numeration is wrong, this corrects it (PDBTMsometimes misses some biomatrix ID numbers)
	bms = sorted(list(PDB_dict['BIOMATRIX'].keys()))
	if bms != list(range(len(bms))):
		for n, bmi in enumerate(bms):
			PDB_dict['BIOMATRIX'][n] = copy.deepcopy(PDB_dict['BIOMATRIX'][bmi])


	# 5. Includes path of the parsed structure
	PDB_dict['PATH'] = new_path
	PDB_dict['HEADER_PATH'] = new_path
	print(struct, PDB_dict['PATH'], 'FROM', pdbfrom, 'BIOLOGICAL_UNIT', PDB_dict['BIOLOGICAL_UNIT'], 'CHAINS', PDB_dict['CHAINS'], 'BIOMATRIX', PDB_dict['BIOMATRIX'], 'Sequence', [(x, PDB_dict[x]['RESNAMES']) for x in PDB_dict['CHAINS']])


	# Is it a full structure? Then, it might have identical sets of chains
	if identical_chain_sets:
		PDB_dict['IDENTICAL_SETS'] = identical_chain_sets
	else:
		PDB_dict['IDENTICAL_SETS'] = {0: list(PDB_dict['CHAINS'])}

	if hkeys:
		return PDB_dict, header_keys
	if header_from:
		for k in header_keys:
			PDB_dict[k] = header_from[k]
	return PDB_dict



### 2. Error logging

# Printing logs
def print_log(message, thr_log_status='ERROR', log_filename=''):
	"""
	Function to print log messages.
	"message" is a 3-tuple (STATUS, FUNC_NAME, TEXT). 
	STATUS must be one of the following:
	  DEBUG for debugging purposes
	  NOTICE for tricky things happening under the hood
	  WARNING if something could create an error
	  ERROR if a process fails
	  CRITICAL if there is a contradiction and the program has to stop
	FUNC_NAME is the name of the function calling this function
	TEXT is a string of text, that can contain any kind of printable
	structure and any number of lines.
	"thr_log_status" determines which stati are also printed in stdout.
	If thr_log_status="ERROR", stati ERROR and CRITICAL are printed.
	"""
	if not log_filename:
		if multiprocessing.current_process().name == 'MainProcess':
			log_filename = 'DEBUG_' + str(os.getpid()) + '_log.txt'
		else:
			log_filename = 'DEBUG_' + str(os.getppid()) + '_log.txt'

	log_stati = ['DEBUG', 'NOTICE', 'WARNING', 'ERROR', 'CRITICAL']

	# Checks over print_log itself and checks fixed arguments
	if thr_log_status not in log_stati:
		print_log(('CRITICAL', print_log.__name__, 'threshold_log_status \'{0}\' is not recognized. List of available log_stati: {1}'.format(thr_log_status, log_stati)))
	if len(message) != 3:
		print_log(('CRITICAL', print_log.__name__, 'message argument should be a 3-tuple: (log_status, function_name, message_text). Here, message was {0}'.format(message)))

	# Defines variables
	log_status = message[0]
	function_logged = message[1]
	message_texts = message[2].split('\n')
	ts = time.time()
	file_time_stamp = datetime.datetime.fromtimestamp(ts).strftime('%Y%m%d_%H%M%S')
	time_stamp = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
	if len(function_logged) > 30:
		formatted_function_logged = function_logged[:27] + '...'
	else:
		formatted_function_logged = function_logged

	# Checks arguments 
	if log_status not in log_stati:
		print_log(('CRITICAL', print_log.__name__, 'message log_status argument \'{0}\' is not recognized. List of available log_stati: {1}'.format(log_status, log_stati)))

	# Formats the output
	thr_ilog = log_stati.index(thr_log_status)
	msg_ilog = log_stati.index(log_status)	
	formatted_info = '{0:20}\t{1:10}\t{2:30}\t'.format(time_stamp, log_status, formatted_function_logged)
	init = True 
	if not log_filename:
		log_filename = './log__' + function_logged + '__' + file_time_stamp + '.txt'
	log_file = open(log_filename, 'a')
	for message_text in message_texts:
		nc = 0
		while nc < len(message_text):
			if len(message_text[nc:]) > 60:
				last_space_i = message_text[nc:nc+60].rfind(' ')
				if last_space_i == -1:
					new_nc = nc+60
				else:
					new_nc = nc+last_space_i
			else:
				new_nc = nc + len(message_text[nc:])
			formatted_text = formatted_info + message_text[nc:new_nc].lstrip() + '\n'
			if init:
				formatted_info = ' '*20 + '\t' + ' '*10 + '\t' + ' '*30 + '\t'
				init = False
	
			log_file.write(formatted_text)
			if msg_ilog >= thr_ilog:
				print(formatted_text[:-1])	# It already has a '\n'
			nc = new_nc


def test_print_log():
	"""
	Test for the log function
	"""
	print_log(('DEBUG', __name__, 'Hi!\n As a test, I will print a small message and a list of numbers passed as an argument. Bye! {0}'.format(list(range(100)))), thr_log_status='DEBUG')


### 3. Downloading

# Downloading FASTAs
def download_FASTAs(fasta_list, out_dir='./', skip_if_present=False, read_not_done_filename='', write_not_done_filename='', thr_log_status='ERROR'):
	"""
	Downloads a list of fasta files from the PDB, and can annotate when
	the files fail the download.
	"""
	if out_dir == './' and len(fasta_list) > 10:
		print_log(('WARNING', download_FASTAs.__name__, 'This function is going to download {0} files in this folder ({1})'.format(len(fasta_list), os.getcwd())))

	not_done = set()
	if read_not_done_filename and os.path.exists(read_not_done_filename):
		not_done_file = open(read_not_done_filename, 'r')
		text = not_done_file.read().split('\n')
		not_done_file.close()

		for line in text:
			if not line:
				continue
			fields = line.split()
			not_done.add((fields[0], fields[1]))

	for pdbname in fasta_list:
		if pdbname in [x[0] for x in not_done]:
			continue
		fasta_filename = out_dir + pdbname + '.fa'
		if skip_if_present and os.path.exists(fasta_filename) and os.path.getsize(fasta_filename) > 0:
			print_log(('NOTICE', download_FASTAs.__name__, 'FASTA file {0} already present in output directory {1}: download skipped (to overwrite, use default option \'skip_if_present=False\')'.format(fasta_filename, out_dir)), thr_log_status)
			continue
		fasta_address = 'https://www.rcsb.org/pdb/download/viewFastaFiles.do?structureIdList={0}\&compressionType=uncompressed'.format(pdbname.upper())
		os.system("wget {0} -O {1} 2>&1 > /dev/null".format(fasta_address, fasta_filename))
		if not os.path.exists(fasta_filename) or os.path.getsize(fasta_filename) == 0:
			print_log(('ERROR', download_FASTAs.__name__, 'FASTA file {0} could not be downloaded from {1}'.format(pdbname, fasta_address)), thr_log_status)
			not_done.add((pdbname, fasta_filename))

	if write_not_done_filename:
		not_done_file = open(write_not_done_filename, 'w')
		for pfam, pfam_filename in not_done:
			not_done_file.write(pfam + '\t' + pfam_filename + '\n')
		not_done_file.close()


# Downloading PDBs
def download_PDBs(pdb_list, out_dir='./', skip_if_present=False, read_not_done_filename='', write_not_done_filename='', thr_log_status='ERROR'):
	"""
	Downloads a list of PDB files from the PDB, and can annotate when
	the files fail the download.
	"""
	if out_dir == './' and len(pdb_list) > 10:
		print_log(('WARNING', download_PDBs.__name__, 'This function is going to download {0} files in this folder ({1})'.format(len(pdb_list), os.getcwd())))

	not_done = set()
	if read_not_done_filename and os.path.exists(read_not_done_filename):
		not_done_file = open(read_not_done_filename, 'r')
		text = not_done_file.read().split('\n')
		not_done_file.close()

		for line in text:
			if not line:
				continue
			fields = line.split()
			not_done.add((fields[0], fields[1]))

	for pdbname in pdb_list:
		if pdbname in [x[0] for x in not_done]:
			continue
		pdb_filename = out_dir + pdbname + '.pdb'
		if skip_if_present and os.path.exists(pdb_filename) and os.path.getsize(pdb_filename) > 0:
			print_log(('NOTICE', download_PDBs.__name__, 'PDB file {0} already present in output directory {1}: download skipped (to overwrite, use default option \'skip_if_present=False\')'.format(pdb_filename, out_dir)), thr_log_status)
			continue
		pdb_address = 'https://files.rcsb.org/view/{0}.pdb'.format(pdbname.upper())
		os.system("wget {0} -O {1} 2>&1 > /dev/null".format(pdb_address, pdb_filename))
		if not os.path.exists(pdb_filename) or os.path.getsize(pdb_filename) == 0:
			print_log(('ERROR', download_PDBs.__name__, 'PDB file {0} could not be downloaded from {1}'.format(pdbname, pdb_address)), thr_log_status)
			not_done.add((pdbname, pdb_filename))

	if write_not_done_filename:
		not_done_file = open(write_not_done_filename, 'w')
		for pfam, pfam_filename in not_done:
			not_done_file.write(pfam + '\t' + pfam_filename + '\n')
		not_done_file.close()


# Downloading Pfam alignments and hmms
def download_Pfams(pfam_list, out_dir='./', skip_if_present=False, read_not_done_filename='', write_not_done_filename='', thr_log_status='ERROR'):
	"""
	Downloads a list of Pfam alignments and HMMs, and can annotate when
	the files fail the download.
	"""
	if out_dir == './' and len(pdb_list) > 10:
		print_log(('WARNING', download_Pfams.__name__, 'This function is going to download {0} files in this folder ({1})'.format(len(pdb_list), os.getcwd())))

	not_done = set()
	if read_not_done_filename and os.path.exists(read_not_done_filename):
		not_done_file = open(read_not_done_filename, 'r')
		text = not_done_file.read().split('\n')
		not_done_file.close()

		for line in text:
			if not line:
				continue
			fields = line.split()
			not_done.add((fields[0], fields[1]))

	for pfam in pfam_list:
		if pfam in [x[0] for x in not_done]:
			continue
		pfam_filename = out_dir + pfam + '.hmm'
		if skip_if_present and os.path.exists(pfam_filename) and os.path.getsize(pfam_filename) > 0:
			print_log(('NOTICE', download_Pfams.__name__, 'Pfam file {0} already present in output directory {1}: download skipped (to overwrite, use default option \'skip_if_present=False\')'.format(pfam_filename, out_dir)), thr_log_status)
		else:
			pfam_address = 'https://pfam.xfam.org/family/{0}/hmm'.format(pfam)
			os.system("wget --no-check-certificate {0} -O {1}".format(pfam_address, pfam_filename))
			if not os.path.exists(pfam_filename) or os.path.getsize(pfam_filename) == 0:
				print_log(('ERROR', download_Pfams.__name__, 'Pfam file {0} could not be downloaded from {1}'.format(pfam_filename, pfam_address)), thr_log_status)
				not_done.add((pfam, pfam_filename))
		pfam_filename = out_dir + pfam + '.sth'
		pfam_fasta_filename = out_dir + pfam + '.fasta'
		if skip_if_present and (os.path.exists(pfam_filename) and os.path.getsize(pfam_filename) > 0 and os.path.exists(pfam_fasta_filename) and os.path.getsize(pfam_fasta_filename) > 0):
			print_log(('NOTICE', download_Pfams.__name__, 'Pfam file {0} already present in output directory {1}: download skipped (to overwrite, use default option \'skip_if_present=False\')'.format(pfam_filename, out_dir)), thr_log_status)
		else:
			if not (os.path.exists(pfam_filename) and os.path.getsize(pfam_filename) > 0):
				pfam_address = 'https://pfam.xfam.org/family/{0}/alignment/full'.format(pfam)
				os.system("wget --no-check-certificate {0} -O {1}".format(pfam_address, pfam_filename))
				if not os.path.exists(pfam_filename) or os.path.getsize(pfam_filename) == 0:
					print_log(('ERROR', download_Pfams.__name__, 'Pfam file {0} could not be downloaded from {1}'.format(pfam_filename, pfam_address)), thr_log_status)
					not_done.add((pfam, pfam_filename))
				else:
					sth2fasta(pfam_filename, pfam_fasta_filename)
			else:
				sth2fasta(pfam_filename, pfam_fasta_filename)
				if not os.path.exists(pfam_fasta_filename) or os.path.getsize(pfam_fasta_filename) == 0:
					print_log(('ERROR', download_Pfams.__name__, 'Pfam stokholm file {0} could not be formatted in fasta'.format(pfam_filename)), thr_log_status)
					not_done.add((pfam, pfam_filename))

	if write_not_done_filename:
		not_done_file = open(write_not_done_filename, 'w')
		for pfam, pfam_filename in not_done:
			not_done_file.write(pfam + '\t' + pfam_filename + '\n')
		not_done_file.close()
	return not_done


### 4. Format editing

# Select one fasta chain
def select_fasta_chain(whole_fasta_filename, chain_fasta_filename, six_pdbcode):
	"""
	From a proper fasta file, prints a file with only one chain.
	"""
	whole_fasta_file = open(whole_fasta_filename, 'r')
	text = whole_fasta_file.read().split('\n')
	whole_fasta_file.close()

	chain_fasta_file = open(chain_fasta_filename, 'w')
	fasta_sequence = ''
	to_copy = False
	for line in text:
		if not line:
			continue
		if line[0] == '>':
			if line[1:5] == six_pdbcode[:4].upper() and line[6] == six_pdbcode[5]:
				to_copy = True
				chain_fasta_file.write(line + '\n')
			else:
				to_copy = False
		else:
			if to_copy:
				chain_fasta_file.write(line + '\n')
				fasta_sequence += line.strip()
	chain_fasta_file.close()
	return fasta_sequence


# Select one PDB chain
def select_pdb_chain(whole_pdb_filename, chain_pdb_filename, six_pdbcode):
	"""
	From a proper PDB file, prints a file with only one chain.
	"""
	whole_pdb_file = open(whole_pdb_filename, 'r')
	text = whole_pdb_file.read().split('\n')
	whole_pdb_file.close()

	chain_pdb_file = open(chain_pdb_filename, 'w')
	fasta_sequence = ''
	to_copy = False
	for line in text:
		if not line:
			continue
		if line.startswith('ATOM') and line[21] == six_pdbcode[5]:
			chain_pdb_file.write(line + '\n')
			# TO DO: retrieve fasta sequence here! From structure, not from SEQRES. Structure is the only thing there must be in a PDB.
	chain_pdb_file.close()
	return fasta_sequence


# Converting Stockholm sequences to fasta format
def sth2fasta(sth_filename, fasta_filename, fmt='utf-8', forget_aln=False):
	"""
	Converts Stockholm files to fasta files.
	The encoding format can be changed using the "fmt" argument.
	The "forget_aln" argument writes the sequence in standard, unaligned fasta.
	"""
	fasta_file = open(fasta_filename, 'w')
	with codecs.open(sth_filename, 'r', fmt) as sth_file:
		for line in sth_file:
			if not line:
				continue
			if line.strip().startswith("#") or line.strip().startswith("//"):
				continue
			fields = line.split()
			fasta_file.write(">" + fields[0] + '\n')
			if forget_aln:
				fasta_file.write(fields[1].replace('-','').replace('.','').upper() + '\n')
			else:
				fasta_file.write(fields[1] + '\n')
	fasta_file.close()


# Writing tables
def write_table(entries, table_filename, tformat='tsv', thr_log_status='ERROR'):
	"""
	From a list of entries (tuples all having the same format) it compiles a tsv
	or csv file (the two options supported by the argument "tformat" at the moment)
	"""
	tformats = {'tsv' : '\t', 'csv' : ','}
	if tformat not in tformats:
		print_log(('ERROR', write_table.__name__, 'the table\'s format {0} is not supported. Supported formats: {1}'.format(tformat, [x for x in tformats])), thr_log_status)

	if not table_filename:
		for e in entries:
			s = ''
			for f in e:
				s += str(f)+'\t'
			print(s)
		return

	table_file = open(table_filename, 'w')
	for e in entries:
		for f in e[:-1]:
			table_file.write("{0}{1}".format(f, tformats[tformat]))
		table_file.write("{0}\n".format(e[-1]))
	table_file.close()
	
	return


# Reading tables
def read_table(table_filename, tformat='tsv', thr_log_status='ERROR'):
	"""
	Reads a tsv or csv-formatted table and returns a list of entries (lists)
	"""
	tformats = {'tsv' : '\t', 'csv' : ','}
	if tformat not in tformats:
		print_log(('ERROR', write_table.__name__, 'the table\'s format {0} is not supported. Supported formats: {1}'.format(tformat, [x for x in tformats])), thr_log_status)

	table_file = open(table_filename, 'r')
	text = table_file.read().split('\n')
	table_file.close()

	entries = []
	for line in text:
		if not line:
			continue
		fields = line.split(tformats[tformat])
		entries.append(fields)
		
	return entries


# Creating an SQL table
def make_sql_table(table_sql_filename, table_entries=[], table_text_filename='', tformat='tsv', thr_log_status='ERROR'):
	"""
	From a list of entries (tuples passed in the argument "table_entries") it
	creates an SQL table.
	Alternatively it can also parse a csv or tsv table in a text file.
	"""
	tformats = {'tsv' : '\t', 'csv' : ','}
	if tformat not in tformats:
		print_log(('ERROR', make_sql_table.__name__, 'in table {0}, the table\'s format {1} is not supported. Supported formats: {2}'.format(table_filename, tformat, [x for x in tformats])), thr_log_status)

	if not table_entries:
		if not table_text_filename or not os.path.exists(table_text_filename):
			print_log(('ERROR', make_sql_table.__name__, 'while building sql table {0} the table entries were not found and the table text file {1} was not found'.format(table_sql_filename, table_text_filename)), thr_log_status)
		else:
			table_file = open(table_text_filename, 'r')
			text = table_file.read().split('\n')
			table_file.close()

			table_entries = []
			declared_fields = text[0].split(tformats[tformat])
			entries += [declared_fields]
			for nl, line in enumerate(text):
				if not line:
					return
				fields = line.split(tformats[tformat])
				if len(fields) != len(declared_fields):
					print_log(('ERROR', make_sql_table.__name__, 'in table {0}, line {1}, there are {2} fields instead of {3} (declared table format: {4})'.format(table_filename, nl, len(fields), len(declared_fields), tformat)), thr_log_status)
				table_entries.append(fields)

	connection = sqlite3.connect(table_sql_filename)
	cursor = connection.cursor()

	sql_types = make_domdom_contact_table('', {}, {}, {}, sql_types=True)
	sql_command = 'CREATE TABLE main ('
	for nf, f in enumerate(table_entries[0]):
		sql_command +=  str(f) + ' ' + sql_types[nf] + ', '
	sql_command = sql_command[:-2] + ');'

	cursor.execute(sql_command)

	for e in table_entries[1:]:
		format_str = 'INSERT INTO main ('
		for df in table_entries[0]:
			format_str += df + ', '
		format_str = format_str[:-2] + ') VALUES ('
		for f in e:
			format_str += '"' + str(f) + '", '
		format_str = format_str[:-2] + ');'
		cursor.execute(format_str)

	connection.commit()
	connection.close()

	return 
