#!/usr/bin/ python

import sys

def getRemainingSites(allSites, selectedSites, symbol):
	for item in selectedSites:
		del allSites[allSites.index(item)]
	
	return allSites, [symbol] * len(allSites)

def getUniqueList(seq):
    return {}.fromkeys(seq).keys()

input_sites = sys.argv[1]

tmpSites = input_sites.split('-')
user_sites = [[int(x) for x in tmpList.split(',')] for tmpList in tmpSites]

Bush = [121, 124, 133, 135, 138, 142, 145, 156, 158, 186, 190, 193, 194, 197, 201, 226, 262, 275]
remaining_sites, remaining_sites_ids = getRemainingSites(range(335), Bush, '_')
remaining_sites.extend(Bush)
remaining_sites_ids.extend(['B'] * len(Bush))
Bush_dict = dict(zip(remaining_sites, remaining_sites_ids))

Fitch = [138, 145, 156, 186, 193, 226]
remaining_sites, remaining_sites_ids = getRemainingSites(range(335), Fitch, '_')
remaining_sites.extend(Fitch)
remaining_sites_ids.extend(['F'] * len(Fitch))
Fitch_dict = dict(zip(remaining_sites, remaining_sites_ids))


Huang2009 = [62, 193, 160, 137, 260, 197, 156, 158, 278, 189, 145]
remaining_sites, remaining_sites_ids = getRemainingSites(range(335), Huang2009, '_')
remaining_sites.extend(Huang2009)
remaining_sites_ids.extend(['H'] * len(Huang2009))
Huang2009_dict = dict(zip(remaining_sites, remaining_sites_ids))

Lee2007_total = [2, 82, 308, 88, 92, 299, 121, 126, 135, 144, 145, 155, 156, 157, 158, 160, 173, 174, 188, 201, 213, 230, 240, 247, 262, 276]
remaining_sites, remaining_sites_ids = getRemainingSites(range(335), Lee2007_total, '_')
remaining_sites.extend(Lee2007_total)
remaining_sites_ids.extend(['L'] * len(Lee2007_total))
Lee2007_total_dict = dict(zip(remaining_sites, remaining_sites_ids))

Lee2007_pos = [82, 92, 121, 135, 144, 145, 155, 156, 157, 158, 160, 173, 174, 188, 240, 247, 276]
remaining_sites, remaining_sites_ids = getRemainingSites(range(335), Lee2007_pos, '_')
remaining_sites.extend(Lee2007_pos)
remaining_sites_ids.extend(['L'] * len(Lee2007_pos))
Lee2007_pos_dict = dict(zip(remaining_sites, remaining_sites_ids))


Liao2008 = [82, 92, 121, 124, 129, 135, 144, 145, 155, 156, 157, 158, 160, 173, 174, 188, 189, 240, 273, 276]
remaining_sites, remaining_sites_ids = getRemainingSites(range(335), Liao2008, '_')
remaining_sites.extend(Liao2008)
remaining_sites_ids.extend(['L'] * len(Liao2008))
Liao2008_dict = dict(zip(remaining_sites, remaining_sites_ids))

Pond2008 = [45, 135, 145, 155, 158, 229, 248, 331]
remaining_sites, remaining_sites_ids = getRemainingSites(range(335), Pond2008, '_')
remaining_sites.extend(Pond2008)
remaining_sites_ids.extend(['P'] * len(Pond2008))
Pond2008_dict = dict(zip(remaining_sites, remaining_sites_ids))


onSurface = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 21, 22, 23, 24, 25, 27, 31, 32, 33, 35, 40, 41, 45, 46, 47, 48, 49, 50, 53, 54, 55, 57, 60, 62, 63, 75, 77, 78, 80, 82, 83, 85, 91, 92, 93, 94, 95, 96, 103, 104, 105, 106, 114, 119, 121, 122, 124, 126, 128, 129, 131, 132, 133, 134, 135, 137, 140, 141, 142, 143, 144, 145, 149, 156, 157, 158, 159, 160, 162, 163, 165, 167, 169, 171, 172, 173, 175, 187, 188, 189, 190, 192, 193, 194, 196, 197, 198, 199, 207, 208, 209, 211, 212, 214, 222, 224, 225, 226, 239, 240, 255, 259, 261, 262, 263, 264, 269, 271, 273, 274, 275, 276, 277, 278, 280, 289, 290, 291, 296, 299, 307, 308, 310, 312, 313, 315, 321, 325, 326, 327, 328 ] 
remaining_sites, remaining_sites_ids = getRemainingSites(range(335), onSurface, '_')
remaining_sites.extend(onSurface)
remaining_sites_ids.extend(['v'] * len(onSurface))
onSurface_dict = dict(zip(remaining_sites, remaining_sites_ids))

# -----------------------------------------------------------------------------
#Hensley, S. E., Hickman, H. D., Jayaraman, A., Viswanathan, K., et al. (2009). Hemagglutinin receptor binding avidity drives influenza A virus antigenic drift.
hensley = [156,158,246] #avidity changing
remaining_sites, remaining_sites_ids = getRemainingSites(range(335), hensley, '_')
remaining_sites.extend(hensley)
remaining_sites_ids.extend(['H'] * len(hensley))

hensley_dict = dict(zip(remaining_sites, remaining_sites_ids))
# -----------------------------------------------------------------------------


# Gamblin, S. J., Haire, L. F., Russell, R. J., Stevens, D. J., Xiao, B., Ha, Y., Vasisht, N., et al. (2004). The structure and receptor binding properties of the 1918 influenza hemagglutinin. Science (New York, N.Y.),
skehel = range(189,200) # 190 helix
skehel.extend(range(133,139)) # 130 loop
skehel.extend(range(224,229)) # 220 loop
skehel.extend([98,153,155,183]) # center
remaining_sites, remaining_sites_ids = getRemainingSites(range(335), skehel, '_')
remaining_sites.extend(skehel)
remaining_sites_ids.extend(['S'] * len(skehel))
skehel_dict = dict(zip(remaining_sites, remaining_sites_ids))

# -----------------------------------------------------------------------------

Xia2009 = [78, 122, 188, 207, 242, 276]
Xia2009.extend([2, 53, 137, 213, 260, 244, 50])
Xia2009.extend([54, 133, 143, 156, 160, 172, 197, 217, 121])
Xia2009.extend([124, 135, 145, 157, 189, 190, 196, 226, 262])
#Xia2009.extend([124, 135, 145, 157, 189, 190, 196, 226, 262, 276])
Xia2009.extend([25, 57, 75, 83, 131, 142, 144, 155, 186, 222, 225, 227, 63, 82, 94, 126, 192, 202, 299])

remaining_sites, remaining_sites_ids = getRemainingSites(range(335), Xia2009, '_')
remaining_sites.extend(Xia2009)
remaining_sites_ids.extend(['X'] * len(Xia2009))
Xia2009_dict = dict(zip(remaining_sites, remaining_sites_ids))

# -----------------------------------------------------------------------------
#Matrosovitch 2000: Early Alterations of the Receptor-Binding Properties of H1, H2, and H3 Avian Influenza Virus Hemagglutinins after Their Introduction into Mammals
matrosovich = range(131,149)
matrosovich.extend(range(152,157))
matrosovich.extend(range(183,196))
matrosovich.extend(range(218,231))
remaining_sites, remaining_sites_ids = getRemainingSites(range(335), matrosovich, '_')
remaining_sites.extend(matrosovich)
remaining_sites_ids.extend(['M'] * len(matrosovich))
matrosovich_dict = dict(zip(remaining_sites, remaining_sites_ids))

# -----------------------------------------------------------------------------
#Aoyama, T., Nobusawa, E., & Kato, H. (1991). Comparison of Complete Amino Acid Sequences among 13 Serotypes of Hemagglutinins and Receptor-Binding Properties of Influenza A Viruses Indirect immunofluorescence.
nobusawa = [98,153,155,183,190,194,195] # Several amino acids known to form the re- ceptor-binding cavity and the second shell in H3 HA (H3 numbering: Tyr-98, Gly-134, Ala-l 38, Trp-153, and His-l 83, Tyr-195) (Weis et al., 1988) are conserved all through the 13 serotypes
nobusawa.extend(range(224,230)) # there are several amino acid residues which compose the receptor-binding site but are not conserved. A vari- ant of H3 HA from which residues HA,224-230 con- sisting the left edge of the receptor-binding site are deleted retained receptor-binding activity (Daniels et a/., 1987). To
nobusawa.extend(range(134,138)) # . While crucial amino acids for receptor binding are conserved, only Gly-134 is conserved on the right edge (134-l 38) of the receptor- binding site all through the 13 serotype HAS

remaining_sites, remaining_sites_ids = getRemainingSites(range(335), nobusawa, '_')
remaining_sites.extend(nobusawa)
remaining_sites_ids.extend(['O'] * len(nobusawa))

nobusawa_dict = dict(zip(remaining_sites, remaining_sites_ids))

# -----------------------------------------------------------------------------
# epitope_sites = [[122, 124, 126, 130, 131, 132, 133, 135, 137, 138, 140, 142, 143, 144, 145, 146, 150, 152, 168],[128, 129, 155, 156, 157, 158, 159, 160, 163, 164, 165, 186, 187, 188, 189, 190, 192, 193, 194, 196, 197, 198],[44, 45, 46, 47, 48, 50, 51, 53, 54, 273, 275, 276, 278, 279, 280, 294, 297, 299, 300, 304, 305, 307, 308, 309, 310, 311, 312], [96, 102, 103, 117, 121, 167, 170, 171, 172, 173, 174, 175, 176, 177, 179, 182, 201, 203, 207, 208, 209, 212, 213, 214, 215, 216, 217, 218, 219, 226, 227, 228, 229, 230, 238, 240, 242, 244, 246, 247, 248], [57, 59, 62, 63, 67, 75, 78, 80, 81, 82, 83, 86, 87, 88, 91, 92, 94, 109, 260, 261, 262, 265]]  by sites
epitope_sites = [122, 124, 126, 130, 131, 132, 133, 135, 137, 138, 140, 142, 143, 144, 145, 146, 150, 152, 168, 128, 129, 155, 156, 157, 158, 159, 160, 163, 164, 165, 186, 187, 188, 189, 190, 192, 193, 194, 196, 197, 198, 44, 45, 46, 47, 48, 50, 51, 53, 54, 273, 275, 276, 278, 279, 280, 294, 297, 299, 300, 304, 305, 307, 308, 309, 310, 311, 312, 96, 102, 103, 117, 121, 167, 170, 171, 172, 173, 174, 175, 176, 177, 179, 182, 201, 203, 207, 208, 209, 212, 213, 214, 215, 216, 217, 218, 219, 226, 227, 228, 229, 230, 238, 240, 242, 244, 246, 247, 248, 57, 59, 62, 63, 67, 75, 78, 80, 81, 82, 83, 86, 87, 88, 91, 92, 94, 109, 260, 261, 262, 265] # full list
epitope_sites_ids = ['A'] * 19
epitope_sites_ids.extend(['B'] * 22)
epitope_sites_ids.extend(['C'] * 27)
epitope_sites_ids.extend(['D'] * 41)
epitope_sites_ids.extend(['E'] * 22)
remaining_sites, remaining_sites_ids = getRemainingSites(range(335), epitope_sites, '_')
remaining_sites.extend(epitope_sites)
remaining_sites_ids.extend(epitope_sites_ids)
epitope_sites_dict = dict(zip(remaining_sites,remaining_sites_ids))

# Neumann 2006: Host range restriction and pathogenicity in the context of influenza pandemic.
# We summarize current knowledge of viral factors that determine host range restriction and pathogenicity of influenza A viruses.
neumann=[226,228]
remaining_sites, remaining_sites_ids = getRemainingSites(range(335), neumann, '_')
remaining_sites.extend(neumann)
remaining_sites_ids.extend(['N'] * len(neumann))

neumann_dict = dict(zip(remaining_sites, remaining_sites_ids))

# -----------------------------------------------------------------------------

#dict_list = [epitope_sites_dict, skehel_dict, matrosovich_dict, nobusawa_dict, hensley_dict, neumann_dict, Xia2009_dict]
dict_list = [epitope_sites_dict, skehel_dict, matrosovich_dict, nobusawa_dict, hensley_dict, neumann_dict, Xia2009_dict, Bush_dict, Fitch_dict, Huang2009_dict, Lee2007_total_dict, Lee2007_pos_dict, Liao2008_dict, Pond2008_dict, onSurface_dict]

# -----------------------------------------------------------------------------

type_transitions = [25, 50, 62, 75, 83, 122, 124, 131, 135, 137, 144, 145, 155, 156, 158, 186, 189, 196, 202, 207, 214, 217, 222, 225, 260, 262, 276, 278]
type_transitions_ids = [[10],[3,10],[9],[10],[10],[1],[5],[10],[6,8],[3],[1],[2,6,8],[1,5,10],[9,10],[3,9],[10],[2,5],[9],[10],[1],[7],[2],[10],[10],[3],[8],[9],[2]]

remaining_sites, remaining_sites_ids = getRemainingSites(range(335), type_transitions, [])
remaining_sites.extend(type_transitions)
remaining_sites_ids.extend(type_transitions_ids)

type_transitions_dict = dict(zip(remaining_sites, remaining_sites_ids))

# -----------------------------------------------------------------------------

cluster_transitions = [25, 50, 53, 54, 62, 75, 82, 83, 122, 124, 131, 133, 137, 143, 144, 145, 146, 147, 155, 156, 158, 160, 164, 172, 174, 188, 189, 190, 193, 196, 197, 201, 202, 207, 213, 217, 222, 225, 230, 244, 260, 262, 276, 278]
cluster_transitions_ids = [10],[3,10],[2,3,4],[4],[4,9],[10],[3,4],[10],[1],[5],[10],[4,7],[2,3],[4],[1],[2,6,7,8],[4],[],[1,5,10],[4,7,9,10],[3,9],[4],[2,3],[4],[2,3],[1],[2,5],[7],[2,3],[9],[4],[2,3],[10],[1],[2,3],[2,4],[10],[10],[3],[4],[3],[7],[9],[2]

remaining_sites, remaining_sites_ids = getRemainingSites(range(335), cluster_transitions, [])
remaining_sites.extend(cluster_transitions)
remaining_sites_ids.extend(cluster_transitions_ids)

cluster_transitions_dict = dict(zip(remaining_sites, remaining_sites_ids))

# -----------------------------------------------------------------------------

# Smith, D. J., Lapedes, A. S., De Jong, J. C., Bestebroer, T. M., Rimmelzwaan, G. F., Osterhaus, A. D. M. E., & Fouchier, R. A. M. (2004). Mapping the antigenic and genetic evolution of influenza virus. Science, 305(5682), 371. Retrieved from http://www.sciencemag.org/cgi/content/abstract/305/5682/371
transitions_names = ['EN72','VI75','TX77','BA79','SI87','BE89','BE92','WU95','SY97','FU02']

for i in range(len(user_sites)):
	user_sites[i].sort()
	sys.stdout.write('%i\t%s\t%i' % (i+1, user_sites[i], len(user_sites[i])))
	for current_dict in dict_list:
		sys.stdout.write('\t')
		for item in user_sites[i]:
			sys.stdout.write(current_dict[item])
	
	current_type_transitions = []
	for item in user_sites[i]:
		current_type_transitions.extend(type_transitions_dict[item])
	
	current_type_transitions = getUniqueList(current_type_transitions)
	current_type_transitions.sort()
	
	sys.stdout.write('\t%i\t' % len(current_type_transitions))
	if len(current_type_transitions) > 0:
		for item in current_type_transitions:
			sys.stdout.write('%s,' % transitions_names[item-1])
	
	current_cluster_transitions = []
	for item in user_sites[i]:
		current_cluster_transitions.extend(cluster_transitions_dict[item])
	
	current_cluster_transitions = getUniqueList(current_cluster_transitions)
	current_cluster_transitions.sort()
	
	sys.stdout.write('\t%i\t' % len(current_cluster_transitions))
	if len(current_cluster_transitions) > 0:
		for item in current_cluster_transitions:
			sys.stdout.write('%s,' % transitions_names[item-1])
	
	sys.stdout.write('\n')





