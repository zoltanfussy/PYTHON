import argparse
import pyKEGG

def KO_to_pathways(infile):
	#returns a dictionary of pathways for each KEGG ortholog
	KO2path = {}
	definitions = {}
	categories = "ABCD"
	unwanted = {"09112", # Not included in regular maps
				"09113", # Global maps only
				"09180", # Brite Hierarchies
				"09150", # Organismal Systems
				"09151", # Immune system
				"04640", # Hematopoietic cell lineage
				"04610", # Complement and coagulation cascades
				"04611", # Platelet activation
				"04613", # Neutrophil extracellular trap formation
				"04620", # Toll-like receptor signaling pathway
				"04624", # Toll and Imd signaling pathway
				"04621", # NOD-like receptor signaling pathway
				"04622", # RIG-I-like receptor signaling pathway
				"04623", # Cytosolic DNA-sensing pathway
				"04625", # C-type lectin receptor signaling pathway
				"04650", # Natural killer cell mediated cytotoxicity
				"04612", # Antigen processing and presentation
				"04660", # T cell receptor signaling pathway
				"04658", # Th1 and Th2 cell differentiation
				"04659", # Th17 cell differentiation
				"04657", # IL-17 signaling pathway
				"04662", # B cell receptor signaling pathway
				"04664", # Fc epsilon RI signaling pathway
				"04666", # Fc gamma R-mediated phagocytosis
				"04670", # Leukocyte transendothelial migration
				"04672", # Intestinal immune network for IgA production
				"04062", # Chemokine signaling pathway
				"09152", # Endocrine system
				"04911", # Insulin secretion
				"04910", # Insulin signaling pathway
				"04922", # Glucagon signaling pathway
				"04923", # Regulation of lipolysis in adipocytes
				"04920", # Adipocytokine signaling pathway
				"03320", # PPAR signaling pathway
				"04929", # GnRH secretion
				"04912", # GnRH signaling pathway
				"04913", # Ovarian steroidogenesis
				"04915", # Estrogen signaling pathway
				"04914", # Progesterone-mediated oocyte maturation
				"04917", # Prolactin signaling pathway
				"04921", # Oxytocin signaling pathway
				"04926", # Relaxin signaling pathway
				"04935", # Growth hormone synthesis, secretion and action
				"04918", # Thyroid hormone synthesis
				"04919", # Thyroid hormone signaling pathway
				"04928", # Parathyroid hormone synthesis, secretion and action
				"04916", # Melanogenesis
				"04924", # Renin secretion
				"04614", # Renin-angiotensin system
				"04925", # Aldosterone synthesis and secretion
				"04927", # Cortisol synthesis and secretion
				"09153", # Circulatory system
				"04260", # Cardiac muscle contraction
				"04261", # Adrenergic signaling in cardiomyocytes
				"04270", # Vascular smooth muscle contraction
				"09154", # Digestive system
				"04970", # Salivary secretion
				"04971", # Gastric acid secretion
				"04972", # Pancreatic secretion
				"04976", # Bile secretion
				"04973", # Carbohydrate digestion and absorption
				"04974", # Protein digestion and absorption
				"04975", # Fat digestion and absorption
				"04979", # Cholesterol metabolism
				"04977", # Vitamin digestion and absorption
				"04978", # Mineral absorption
				"09155", # Excretory system
				"04962", # Vasopressin-regulated water reabsorption
				"04960", # Aldosterone-regulated sodium reabsorption
				"04961", # Endocrine and other factor-regulated calcium reabsorption
				"04964", # Proximal tubule bicarbonate reclamation
				"04966", # Collecting duct acid secretion
				"09156", # Nervous system
				"04724", # Glutamatergic synapse
				"04727", # GABAergic synapse
				"04725", # Cholinergic synapse
				"04728", # Dopaminergic synapse
				"04726", # Serotonergic synapse
				"04720", # Long-term potentiation
				"04730", # Long-term depression
				"04723", # Retrograde endocannabinoid signaling
				"04721", # Synaptic vesicle cycle
				"04722", # Neurotrophin signaling pathway
				"09157", # Sensory system
				"04744", # Phototransduction
				"04745", # Phototransduction - fly
				"04740", # Olfactory transduction
				"04742", # Taste transduction
				"04750", # Inflammatory mediator regulation of TRP channels
				"09158", # Development and regeneration
				"04320", # Dorso-ventral axis formation
				"04360", # Axon guidance
				"04361", # Axon regeneration
				"04380", # Osteoclast differentiation
				"09149", # Aging
				"04211", # Longevity regulating pathway
				"04212", # Longevity regulating pathway - worm
				"04213", # Longevity regulating pathway - multiple species
				"09159", # Environmental adaptation
				"04710", # Circadian rhythm
				"04713", # Circadian entrainment
				"04711", # Circadian rhythm - fly
				"04712", # Circadian rhythm - plant
				"04714", # Thermogenesis
				"04626", # Plant-pathogen interaction
				"09160", # Human Diseases
				"09161", # Cancer: overview
				"05200", # Pathways in cancer
				"05202", # Transcriptional misregulation in cancer
				"05206", # MicroRNAs in cancer
				"05205", # Proteoglycans in cancer
				"05204", # Chemical carcinogenesis - DNA adducts
				"05207", # Chemical carcinogenesis - receptor activation
				"05208", # Chemical carcinogenesis - reactive oxygen species
				"05203", # Viral carcinogenesis
				"05230", # Central carbon metabolism in cancer
				"05231", # Choline metabolism in cancer
				"05235", # PD-L1 expression and PD-1 checkpoint pathway in cancer
				"09162", # Cancer: specific types
				"05210", # Colorectal cancer
				"05212", # Pancreatic cancer
				"05225", # Hepatocellular carcinoma
				"05226", # Gastric cancer
				"05214", # Glioma
				"05216", # Thyroid cancer
				"05221", # Acute myeloid leukemia
				"05220", # Chronic myeloid leukemia
				"05217", # Basal cell carcinoma
				"05218", # Melanoma
				"05211", # Renal cell carcinoma
				"05219", # Bladder cancer
				"05215", # Prostate cancer
				"05213", # Endometrial cancer
				"05224", # Breast cancer
				"05222", # Small cell lung cancer
				"05223", # Non-small cell lung cancer
				"09172", # Infectious disease: viral
				"05166", # Human T-cell leukemia virus 1 infection
				"05170", # Human immunodeficiency virus 1 infection
				"05161", # Hepatitis B
				"05160", # Hepatitis C
				"05171", # Coronavirus disease - COVID-19
				"05164", # Influenza A
				"05162", # Measles
				"05168", # Herpes simplex virus 1 infection
				"05163", # Human cytomegalovirus infection
				"05167", # Kaposi sarcoma-associated herpesvirus infection
				"05169", # Epstein-Barr virus infection
				"05165", # Human papillomavirus infection
				"09171", # Infectious disease: bacterial
				"05110", # Vibrio cholerae infection
				"05120", # Epithelial cell signaling in Helicobacter pylori infection
				"05130", # Pathogenic Escherichia coli infection
				"05132", # Salmonella infection
				"05131", # Shigellosis
				"05135", # Yersinia infection
				"05133", # Pertussis
				"05134", # Legionellosis
				"05150", # Staphylococcus aureus infection
				"05152", # Tuberculosis
				"05100", # Bacterial invasion of epithelial cells
				"09174", # Infectious disease: parasitic
				"05146", # Amoebiasis
				"05144", # Malaria
				"05145", # Toxoplasmosis
				"05140", # Leishmaniasis
				"05142", # Chagas disease
				"05143", # African trypanosomiasis
				"09163", # Immune disease
				"05310", # Asthma
				"05322", # Systemic lupus erythematosus
				"05323", # Rheumatoid arthritis
				"05320", # Autoimmune thyroid disease
				"05321", # Inflammatory bowel disease
				"05330", # Allograft rejection
				"05332", # Graft-versus-host disease
				"05340", # Primary immunodeficiency
				"09164", # Neurodegenerative disease
				"05010", # Alzheimer disease
				"05012", # Parkinson disease
				"05014", # Amyotrophic lateral sclerosis
				"05016", # Huntington disease
				"05017", # Spinocerebellar ataxia
				"05020", # Prion disease
				"05022", # Pathways of neurodegeneration - multiple diseases
				"09165", # Substance dependence
				"05030", # Cocaine addiction
				"05031", # Amphetamine addiction
				"05032", # Morphine addiction
				"05033", # Nicotine addiction
				"05034", # Alcoholism
				"09166", # Cardiovascular disease
				"05417", # Lipid and atherosclerosis
				"05418", # Fluid shear stress and atherosclerosis
				"05410", # Hypertrophic cardiomyopathy
				"05412", # Arrhythmogenic right ventricular cardiomyopathy
				"05414", # Dilated cardiomyopathy
				"05415", # Diabetic cardiomyopathy
				"05416", # Viral myocarditis
				"09167", # Endocrine and metabolic disease
				"04930", # Type II diabetes mellitus
				"04940", # Type I diabetes mellitus
				"04950", # Maturity onset diabetes of the young
				"04936", # Alcoholic liver disease
				"04932", # Non-alcoholic fatty liver disease
				"04931", # Insulin resistance
				"04933", # AGE-RAGE signaling pathway in diabetic complications
				"04934", # Cushing syndrome
				"04013", # MAPK signaling pathway - fly
				"04016", # MAPK signaling pathway - plant
				"04011", # MAPK signaling pathway - yeast
				"04341", # Hedgehog signaling pathway - fly
				"04391", # Hippo signaling pathway - fly
				"04392", # Hippo signaling pathway - multiple species
				"04140", # Autophagy - animal
				"04138", # Autophagy - yeast
				"04136", # Autophagy - other
				"04137", # Mitophagy - animal
				"04139", # Mitophagy - yeast
				"04111", # Cell cycle - yeast
				"04112", # Cell cycle - Caulobacter
				"04113", # Meiosis - yeast
				"04214", # Apoptosis - fly
				"05111", # Biofilm formation - Vibrio cholerae
				"02025", # Biofilm formation - Pseudomonas aeruginosa
				"02026", # Biofilm formation - Escherichia coli
				"04550", # Signaling pathways regulating pluripotency of stem cells
				"99994", # Others
				"99998", # Others
				"99999", # Others
				"09194", # Poorly characterized
				"99996", # General function prediction only
				"99997", # Function unknown

				}

	with open(infile) as f:
		for l in f:
			if len(l) == 2:
				continue
			if l[0] not in categories:
				#print(category)
				continue
			category = l[0]
			l = l[1:].strip()
			try:
				if category == "A":
					a_lvl = l.split()[0]
					definitions[a_lvl] = " ".join(l.split()[1:])
				elif category == "B":
					b_lvl = l.split()[0]
					definitions[b_lvl] = " ".join(l.split()[1:])
				elif category == "C":
					c_lvl = l.split()[0]		
					definitions[c_lvl] = " ".join(l.split()[1:])	
				elif category == "D":
					#this is the KEGG ortholog level
					KO = l.split()[0]
					#definition = " ".join(l.split()[1:])
					definitions[KO] = " ".join(l.split()[1:])
					if KO not in KO2path:
						KO2path[KO] = {"A": set(), "B": set(), "C": set()}
					if a_lvl not in unwanted:
						KO2path[KO]["A"].add(a_lvl)
					if b_lvl not in unwanted:
						KO2path[KO]["B"].add(b_lvl)
					if c_lvl not in unwanted:
						KO2path[KO]["C"].add(c_lvl)
			except IndexError:
				print(l)
				quit()
	return definitions, KO2path


def write_table(infile, outfile, definitions, KO2path):
	with open(infile, "rt") as f,\
		 open(outfile, "wt") as result:
		result.write("GeneID\tKEGG-Ortholog\tPATHWAYS\tlevel\tKO_definition\tPWY_definition\n")
		for l in f:
			l = l.split()
			if len(l) == 2:
				geneid = l[0]
				keggid = l[1]
				ko_defline = definitions[keggid]
				pathwaylist_a = list(KO2path[keggid]["A"])
				pathways = ";".join(pathwaylist_a)
				pwy_defline = ";".join([definitions.get(i, "N/A") for i in pathwaylist_a])
				writestring = "{}\t{}\t{}\t{}\t{}\t{}\n".format(geneid, keggid, pathways, "A", ko_defline, pwy_defline)
				result.write(writestring)

				pathwaylist_b = list(KO2path[keggid]["B"])
				pathways = ";".join(pathwaylist_b)
				pwy_defline = ";".join([definitions.get(i, "N/A") for i in pathwaylist_b])
				writestring = "{}\t{}\t{}\t{}\t{}\t{}\n".format(geneid, keggid, pathways, "B", ko_defline, pwy_defline)
				result.write(writestring)

				pathwaylist_c = list(KO2path[keggid]["C"])
				pathways = ";".join(pathwaylist_c)
				pwy_defline = ";".join([definitions.get(i, "N/A") for i in pathwaylist_c])
				writestring = "{}\t{}\t{}\t{}\t{}\t{}\n".format(geneid, keggid, pathways, "C", ko_defline, pwy_defline)
				result.write(writestring)
				#alternatively, write each pathway on a new line, i.e. multiple lines for one geneid
			else:
				geneid = l[0]
				writestring = "{}\t\t\t\t\t\n".format(geneid)
				result.write(writestring)


def main():
	definitions, KO2path = KO_to_pathways("ko00001.keg")

	parser = argparse.ArgumentParser(description='How to use argparse')
	parser.add_argument('-i', '--infile', help='KO file to analyze', required=True)
	#parser.add_argument('-i', '--infile', help='KO file to analyze', default="")
	parser.add_argument('-o', '--outfile', help='KO file to analyze', default="")

	args = parser.parse_args()
	infile = args.infile
	#infile = "/Users/morpholino/OwnCloud/Robolab/KEGG_path/Nfow.ko.txt"
	if args.outfile == "":
		outfile = infile.replace(".ko.txt", "_annot.tsv")
	else:
		outfile = args.outfile

	write_table(infile, outfile, definitions, KO2path)


if __name__ == '__main__':
	main()