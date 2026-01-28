from sys import argv
import pandas as pd

genenames=[]

with open(argv[1]) as f1:
	for f in f1:		
		splitfirst=f.strip().split()
		if splitfirst[0] not in genenames:
			genenames.append(splitfirst[0])
	for i in genenames:
		currentgenelist=[]
		f1.seek(0)
		for f in f1:
			splitagain=f.strip().split()
			if i==splitagain[0]:
				currentgenelist.append(splitagain)
		df=pd.DataFrame(currentgenelist,columns=["source","target","value"])
		df.sort_values(by="value",ascending=False)
		topten=int(len(currentgenelist))//10
		secondcounter=0 
		while secondcounter <= topten:
			printthis=df.loc[secondcounter,:].values.tolist()
			print("	".join(printthis))
			secondcounter+=1
