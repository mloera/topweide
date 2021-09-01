import sys
import re

rc_file = sys.argv[1]

with open( rc_file, 'r' ) as infile:

	txt = infile.readlines()
	print(txt[0].strip() )
	for r in txt[1:]:
	
		r_clean = r.strip().split('\t')

		## Parsing the allele frequencies into floats
		t = [str(int(n.split('/')[0])/int(n.split('/')[1]))  for n in r_clean if re.search('/', n) and not re.search('^[ATGC]', n)]

		## Getting the number of samples (ns)
		ns = int(len(t)/2)

		## Getting the overall major allele for the SNP

		maa = r_clean[4].split('/')[0]

		## Getting the major allele for each sample

		maa_s = [a for a in r_clean[7] ] 



		## Getting the read count per sample and SNP
		t_high = [int(n.split('/')[1])  for n in r_clean if re.search('/', n) and not re.search('^[ATGC]', n)]

		mia_corr = []
		if min(t_high) > int(sys.argv[2]):
			if int(r_clean[3].strip()) == 2:

				## Getting the corrected minor allele frequencies based on the allele state
				for i in range(0, ns):

					if maa_s[i] == maa:
						mia_corr.append(str(1 - float(t[i])))

					else:
						mia_corr.append(str(float(t[i])))

				## Filtering for minimum allele frequency and printing SNP info and min. allele freqs.
				if min([ float(mia) for mia in mia_corr ]) >= float(sys.argv[3]):		
					print('%s\t%s' % ( '\t'.join(r_clean[0:9]), '\t'.join(mia_corr) ) )
