def tss_region(test):
    """take a genes dataframe, calcualte the tss region of each gene
    returns a new dataframe with tss regions and genes name"""
    tss=[]
    for i in range(test.shape[0]):
        if test.iloc[i]['strand']=='-':
            if i+1<test.shape[0] and test.iloc[i+1]['strand']=='-'and test.iloc[i]['seq_ID']==test.iloc[i+1]['seq_ID']: ## both genes are on the '-' strand
                start=test.iloc[i]['end']
                end=test.iloc[i+1]['start']
                tss.append((test.iloc[i]['seq_ID'],start,end, test.iloc[i]['strand'],test.iloc[i]['attributes'].split(';')[0]))
            else: # the next one is on the '+' strand, just go up 100 bases for TSS
                start=test.iloc[i]['end']
                end=start+100
                tss.append((test.iloc[i]['seq_ID'],start,end, test.iloc[i]['strand'],test.iloc[i]['attributes'].split(';')[0]))
        else:## the ones on the '+' strand
            if i>0 and test.iloc[i-1]['strand']=='+'and test.iloc[i]['seq_ID']==test.iloc[i-1]['seq_ID']: ## both on the '+' strand
                start=test.iloc[i-1]['end']
                end=test.iloc[i]['start']
                tss.append((test.iloc[i]['seq_ID'],start,end, test.iloc[i]['strand'],test.iloc[i]['attributes'].split(';')[0]))
            else:
                start=test.iloc[i]['start']-100
                end=test.iloc[i]['start']
                tss.append((test.iloc[i]['seq_ID'],start,end, test.iloc[i]['strand'],test.iloc[i]['attributes'].split(';')[0]))
    test_tss=pd.DataFrame.from_records(tss)
    test_tss.columns=['seqID','tss_start','tss_end','strand','gene_id']
    test_tss['tss_l']=test_tss.tss_end-test_tss.tss_start
    return test_tss
