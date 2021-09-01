import sys
import re

rc_file = sys.argv[1]

with open( rc_file, 'r' ) as infile:

	txt = infile.readlines()
	print(txt[0].strip() )
	for r in txt[1:]:
	
		r_clean = r.strip().split('\t')

		t = [str(int(n.split('/')[0])/int(n.split('/')[1]))  for n in r_clean if re.search('/', n) and not re.search('^[ATGC]', n)]
		t_high = [int(n.split('/')[1])  for n in r_clean if re.search('/', n) and not re.search('^[ATGC]', n)]

		if min(t_high) > int(sys.argv[2]):
				
			print( '%s\t%s' % ( '\t'.join(r_clean[0:9]), '\t'.join(t) ) )
