#!/usr/bin/env python2
# import pymisca.header as pyheader
import pysam
import htseq_ext.htseq_extra as htseq_extra
import HTSeq


import pymisca.vis_util as pyvis
plt = pyvis.plt
import pymisca.ext as pyext
import collections
import numpy as np

def read_gff(gff_filename,
    feature_type=['CDS'],
    id_attribute='Parent',
    additional_attributes=[],
    quiet=0,
    head = -1,
    stranded='yes',
             
            ):
    '''Adapter from HTSeq-count
    '''
    
    if isinstance(feature_type,basestring):
        feature_type = [feature_type]

#     gff_filename = FNAME
    features = HTSeq.GenomicArrayOfSets("auto", stranded != "no")
    gff = HTSeq.GFF_Reader(gff_filename)
    counts = {}    
    attributes = {}
    i = 0
    ids = collections.OrderedDict()
    try:
        for f in gff:
            if f.type in feature_type:
                try:
                    feature_id = f.attr[id_attribute]
                except KeyError:
                    raise ValueError(
                            "Feature %s does not contain a '%s' attribute" %
                            (f.name, id_attribute))
                if stranded != "no" and f.iv.strand == ".":
                    raise ValueError(
                            "Feature %s at %s does not have strand information but you are "
                            "running htseq-count in stranded mode. Use '--stranded=no'." %
                            (f.name, f.iv))
                features[f.iv] += feature_id
                counts[f.attr[id_attribute]] = 0
                attributes[f.attr[id_attribute]] = [
                        f.attr[attr] if attr in f.attr else ''
                        for attr in additional_attributes]
                ids.setdefault( feature_id , [])
                ids[feature_id] += [f.iv]
            i += 1
            if i % 100000 == 0 and not quiet:
                sys.stderr.write("%d GFF lines processed.\n" % i)
                sys.stderr.flush()
            if head >= 0:
                if len(ids) ==head:
                    break
    except:
        sys.stderr.write(
            "Error occured when processing GFF file (%s):\n" %
            gff.get_line_number_string())
        raise
    return features,ids

def dump__ivdqs(ivdq_transcript, BAM):
    if not isinstance(BAM,pysam.AlignmentFile):
        BAM = pysam.AlignmentFile(BAM)

    iv = ivdq_transcript.iv
    reads = BAM.fetch( 
        start=iv.start,
        end=iv.end,
        reference=iv.chrom,)

    ## for debug
    reads_old  = reads= list(reads)    
    ### [FRAG] 
    #### filter for strandness
    reads = (aln for aln in reads if aln.is_reverse == (ivdq_transcript[0].strand=='-') )

    ivdqs_aln = map(pyRiboSeq.htseq_extra.GenomicIntervalDeque.fromAlignment,
                    reads)

    isComp = map(ivdq_transcript.compatibleWithRead, ivdqs_aln)

    print 'Strand of feature: %s' % (ivdq_transcript[0].strand)
    print 'N_aligned:%d'%len(ivdqs_aln)
    print 'Compatible prob:%.4f'% np.mean(isComp)
    return ivdqs_aln

def plotter(ivdq_transcript, ivdqs_aln,ax=None, nMax=1000,colorKey=None,debug = 0,
            upStreamPad=100):
    
    if colorKey is None:
        colorKey = lambda m,r: pyRiboSeq.htseq_comp.model__read__countMismatch(m,r) != 0

    if ax is None:
        figs,axs= plt.subplots(1,1,figsize=[12,6])
        ax = axs
        
#     nMax =1000
#     cmap = plt.get_cmap()
    # ivdqs = [ivdq_transcript,] + ivdqs_aln[:100]
    cmap = pyvis.cmap__defaultCycle()
    ivdqs = [ivdq_transcript,] + ivdqs_aln
    
    if nMax >= 1:
        ivdqs = ivdqs[:nMax] 
    N  = len(ivdqs)

    ivdq_all =  pyRiboSeq.htseq_extra.GenomicIntervalDeque([])
    map(ivdq_all.extend,ivdqs)

    ivdq_all = ivdqs[0]
    
    for i,ivdq in enumerate(ivdqs):
#         arr =
        arr = HTSeq.StepVector.StepVector.create(start_index = ivdq_all.iv.start - upStreamPad)
#         arr = HTSeq.GenomicArray(chroms='auto')
        for x in ivdq:
            if not x.length:
                print (ivdq.cigarstring)
                continue

            arr[x.start:x.end] += 1

#         v_it = arr[ivdq_all.iv].values()
        v_it = arr[ivdq_all.iv.start:ivdq_all.iv.end]
        vals = np.fromiter(
                v_it,
                dtype='float',
            )
        vals[vals==0]= np.nan
        if debug:
            print ('feat_idx',i,'mismatch',
                   pyRiboSeq.htseq_comp.model__read__countMismatch(ivdq_transcript, ivdq , ivdq.iv)
                  )
    
        ax.plot(
            i  + vals,color=cmap(colorKey(ivdq_transcript, ivdq  )),
        )
        
#         color = 1 - ivdq_transcript.compatibleWithRead(ivdq)
#         color = colorKey(  )
#         print (colorKey(ivdq_transcript, ivdq  ))
#         if i == 2:
#             break
    # iv
    #     print np.nanmin(vals,)
    ax.set_ylim(0,None)
    ax.grid(1)
    ax.set_xlabel('Stranded Genomic coordinate (bp)')
    ax.set_ylabel('rank of feature')
#     return ivdqs_aln, ax
    
    return  ax,ivdqs_aln




