import pymisca.tree
# from HTSeq
from HTSeq import GenomicInterval
### lookup table is the ultimate representation of any representable functions
strandDict = pymisca.tree.TreeDict.from_flatPathDict(
        {
        '+/+':'+',
        '+/-':'-',
        '+/.':'.',
        '-/+':'-',
        '-/-':'+',
        '-/.':'.',
        './+':'.',
        './-':'.',
        './.':'.',
        },sep='/')


def model__read__compatibleSpliced(ivdq_model, 
                                   ivdq_read, 
                                   ignoreFalse=True):
    '''
    Usage:
        Decide whether two features overlap significantly
        
    params:
    
    [TODO]:
        Compatibility should be adjusted as CIGAR string in the future
    '''
    ### +++ helper
    def checkCompatible((i,L), iv1,iv2):
        ignoreStart = False
        ignoreEnd = False
        ### only ignore if at 5utr or 3utr
        if (iv1 is ivdq1[0]):# or (iv2 is ivdq2[0]):
            ignoreStart = True
        if (iv1 is ivdq1[-1]):# or (iv2 is ivdq2[-1]):
            ignoreEnd = True

        left = iv2.start >= iv1.start
        right = iv2.end <= iv1.end
        res = True & ( ignoreStart | left) & (ignoreEnd | right )
        assert ignoreFalse or res
        return res

#         left = iv2.start_as_pos >= iv1.start_as_pos
#         rght = iv2.end_as_pos <= iv1.end_as_pos

#         res = (iv2.is_contained_in(iv1)) or (igno)
#         res = True & (ignoreStart or (iv1.start_as_pos == iv2.start_as_pos)) \
#             & (ignoreEnd or (iv1.end_as_pos == iv2.end_as_pos) )

    def findOverlapAndCheck(comp,
                            ivdq1_LeftIter,
                            (i,iv2)
                           ):
        pos = iv2.start_as_pos
        it = (d for d in ivdq1_LeftIter
                    if iv2.overlaps(d['iv']) or d['iv'] is None)

        iv1 = next(it)['iv']

        if iv1 is not None:
#             comp = 0
            if comp is None:
                comp = True
            comp &= checkCompatible((i,L), iv1,iv2)                
        else:
            if iv2 is ivdq2[0]:
                comp = False
#                     assert 0,'Contradict with the overlap() check'
            else:
                #### ivdq2 is longer than ivdq1
#                     comp = False
                assert comp is not None

        return comp,iv1
    ### --- helpers

    #### main function
#     for x in [ivdq_model,ivdq_read]:
#         assert isinstance(x, pyRiboSeq.htseq_extra.GenomicIntervalDeque )
    ivdq1 = ivdq_model
    ivdq2 = ivdq_read
#     if len(ivdq2) > len(ivdq1):
#         ivdq1, ivdq2 = ivdq2, ivdq1
        #### make sure idvq2 is the shorter one    
    
    if not ivdq2.iv.overlaps(ivdq1.iv):        
        return False
    
    
    L = len(ivdq2)
    comp = None

    ivdq1_LeftIter = ivdq1.leftIter(enum=False,copy=True)

    for i, iv2 in enumerate(ivdq2):
        comp,iv1 = findOverlapAndCheck( comp, ivdq1_LeftIter, (i,iv2))  
        if not comp or not iv1:
            break                
    return comp

                
def model__read__countMismatch(ivdq_model, ivdq_aln, 
                               ignoreFalse=True,
                               ignoreModelStart=False,
                               ignoreModelEnd = False,
                               maxMiss = 0,
#                                maxMismatch = 0,
                              ):
    '''
    Usage:
        Count nMiss=(mismatches in bp) between a 
            ivdq_model (usually transcript), and  
            ivdq_aln  (usually aligned short read)
        
    params:
        maxMismatch: aligner will give up if nMiss exceeded this value
        ignoreModelStart: ignore mismatch at the non-stranded start of the model
        ignoreModelEnd: ignore mismatch at the non-stranded end of the model
        
    Note:
        The reported nMiss is NOT clipped by maxMiss.
    [TODO]:
        Compatibility should be adjusted as CIGAR string in the future
    '''
#     for x in [ivdq_model,ivdq_read]:
#         assert isinstance(x, pyRiboSeq.htseq_extra.GenomicIntervalDeque )
    ivdq1 = ivdq_model
    ivdq2 = ivdq_aln
#     if len(ivdq2) > len(ivdq1):
#         ivdq1, ivdq2 = ivdq2, ivdq1
        #### make sure idvq2 is the shorter one    
    if not ivdq2.iv.overlaps(ivdq1.iv):
        return False


#     def checkCompatible((i,L), iv1,iv2):
    def getOverhangLength((i,L), iv1,iv2):
        
        _ignoreStart = ignoreModelStart & (iv1 is ivdq1[0])
        _ignoreEnd   = ignoreModelEnd & (iv1 is ivdq1[-1])

        leftDiff = iv2.start - iv1.start 
        rightDiff = iv2.end - iv1.end
        
        res = ( 1- _ignoreStart ) * max(-leftDiff, 0 ) \
            + ( 1- _ignoreEnd) * max(rightDiff,0)
        
        assert ignoreFalse or res == 0
        return res
    
#         ignoreStart = False
#         ignoreEnd = False
        ### only ignore if at 5utr or 3utr
#         if (iv1 is ivdq1[0]):# or (iv2 is ivdq2[0]):
#             _ignoreStart = True
#         if (iv1 is ivdq1[-1]):# or (iv2 is ivdq2[-1]):
#             _ignoreEnd = True

#         right = iv2.end <= iv1.end
#         leftDiff = 
#         return 
#         res = True \
#             & ( ignoreStart | leftDiff <= 0 ) \
#             & ( ignoreEnd | rightDiff >= 0 )



    def findOverlapAndCheck(nMiss,
                            ivdq1_LeftIter,
                            (i,iv2)
                           ):
#         pos = iv2.start_as_pos
        
        it = (d for d in ivdq1_LeftIter
                    if iv2.overlaps(d['iv']) or d['iv'] is None)

        iv1 = next(it)['iv']

        if iv1 is not None:
#             nMiss += 
            incr = getOverhangLength((i,L), iv1,iv2)                
        else:
#             nMiss += 
            incr = iv2.length
    
#         print(iv2,iv2.length, incr)
        nMiss += incr
        
            #### this segment does not overlap with any of ivdq1
            
#             if iv2 is ivdq2[0]:

#             else:
#                 #### ivdq2 is longer than ivdq1
#                 assert comp is not None

        return nMiss,iv1

    
    L = len(ivdq2)
    nMiss = 0

    ### use the un-stranded iterator for simplicity
    ### also agrees with CIGAR orientation.
    ivdq1_LeftIter = ivdq1.leftIter(enum=False,copy=True)

    for i, iv2 in enumerate(ivdq2):
        nMiss,iv1 = findOverlapAndCheck( nMiss, ivdq1_LeftIter, (i,iv2))  
        if nMiss > maxMiss or not iv1:
            break
            
    return nMiss


def iv__isReverse(iv=None,strand=None):
    if strand is None:
        strand = iv.strand
    return strand =='-'

def iv__strandCompatible(iv1,iv2):
    ### [FRAGILE]
    res = (iv1.strand == iv2.strand)
    return res

def iv__getRelativeStrand(iv1,iv2):
    return strandDict[iv1.strand][iv2.strand]
    
def iv__relativeTo(iv, ivdq_model, ivdq_model_iter, d_model,debug=False):
    '''
    usage:
        helper function that requires caution
    params:
        iv: the GenomicInterval() to be transformed
        d_model: the most element from ivdq_model_iter
        ivdq_model_iter: a stranded iterator usually being ivdq_model.getStrandedIter()
    '''
    iv0 = iv
    iv = iv.copy()
    d = d_model
    
    if 1:
    #### def getReferenceIv(iv, ivdq_model_iter, d):
        while True:            
            if d['iv'] is None:
                #### The iv does not overlap with the model anywhere
                break
            if d['iv'].overlaps(iv):
                break
            d = next(ivdq_model_iter)

    if d['iv'] is None:
        return None

    ##### +++ Core computation +++
    
    ## relative to stranded_start
    iv.start+=  - d['iv'].start_d
    iv.end  +=  - d['iv'].start_d
    
    ## flip negative coords
    if iv__isReverse(iv):
        iv.start, iv.end = iv.end-1, iv.start+1
        iv.start, iv.end = -iv.start, -iv.end
        
    ## shift by skipLen
    iv.start +=  d['skipLen']
    iv.end   +=  d['skipLen']
    
    ##### --- Core computation ---
    
    ### update
    iv.strand = strandDict[iv.strand][ivdq_model.strand]
    iv.chrom = ivdq_model.name
    
    if debug:
        print('[iv__relativeTo]',iv,d)
    
    return iv
    
def model__read__rebase(ivdq_model,
                        ivdq_read,
                        copy    = True,
                        debug   = False,
                        maxMiss = -1):
    '''
    Calculate the coordinate of an ivdq_read with ivdq_model as the reference
    If transforming self, we should get a contiguous region.
    model__read__rebase(ivdq_model, ivdq_model.iv) should flip the strand if reversed
    model__read__rebase(ivdq_model, ivdq_model) should return a contiguous region
    params:
        maxMiss: used to exclude reads. Deactivated by default = -1
        ivdq_model: reference ivdq
        ivdq_read:  ivdq to be rebased
    '''
    assert ivdq_model.strand in '+-.'
    assert iv__strandCompatible(
        ivdq_read, ivdq_model
    ),'Refuses to calculate relative coordinate for non-compatible strand: %s\n%s'%(ivdq_read, ivdq_model)
    ### [TODO]: support opposite-stranded ivdq_model and ivdq_read
    
    if maxMiss >=0:
        nMiss = ivdq_model.countMissWithRead(ivdq_read,
                                             maxMiss=maxMiss)
        if nMiss >=maxMiss:
            print('[WARN] read is discarded due to mismatch:%s'%(ivdq_read.iv))
            return None
    
    if copy:
        ivdq_read0 = ivdq_read
        ivdq_read  = ivdq_read.copy()
        ivdq_read.last = ivdq_read0
        ### .last will not be set unless copied
        
    ivdq_read_iter  = ivdq_read.getStrandedIter()
    ivdq_model_iter = ivdq_model.getStrandedIter() 
    
    
    ivdq_read.clear()
    if iv__isReverse(ivdq_read):
        ivdq_read.strand = '+'
        
    d_model  = next(ivdq_model_iter)
    for d_read in ivdq_read_iter:
        if d_read['iv'] is None:
            break
        res = iv__relativeTo( d_read['iv'], ivdq_model, ivdq_model_iter, d_model,debug=debug)
        if res is None:
            ### ran out of ivdq_model_iter already. break
            break
        else:
#         if res is not None:
            #### assume output is always forward stranded
            ivdq_read.append(res)
                   
    return ivdq_read



def iv__null():
    iv = GenomicInterval('NA',0,0,strand='.')
    return iv
def iv__intersect(iv1,iv2):
    if not iv1.overlaps(iv2):
        return iv__null()
    else:
        iv = GenomicInterval(
            chrom = iv1.chrom,
            strand = iv1.strand,
            start=max(iv1.start, iv2.start),
            end = min(iv1.end, iv2.end)
        )
        return iv
    
    
def iv__fastaDict__getSequence(iv,fastaDict,seq=None,castFunc=None):
    if castFunc is None:
        castFunc = lambda x:x.upper()
    if seq is None:
        seq = fastaDict[iv.chrom]
    seq = seq[iv.start:iv.end]
    seq = castFunc(seq)
    if iv.strand in ['+','.']:
        pass
    elif iv.strand in ['-']:
#         seq = seq[::-1]
        seq = seq__reverseComplement(seq)
    return seq


def ivdq__fastaDict__getSequence(ivdq, fastaDict,seq=None,castFunc=None):
    if seq is None:
        seq = fastaDict[ivdq[0].chrom]
#     ivdq.sortBy(lambda x:x.start * (x.strand=='-'))
    ivdq = list(ivdq)
    ivdq = sorted(ivdq, key = lambda x:x.start * [1,-1][(x.strand=='-')])

    seqs = [iv__fastaDict__getSequence(iv,fastaDict, seq=seq, castFunc=castFunc) for iv in ivdq]
    return ''.join(seqs)


DNA_COMPLEMENTER = {'A':'T','C':'G','T':'A','G':'C','N':'N'}
def seq__complement(seq,mapper= None):
    if mapper is None:
        mapper = DNA_COMPLEMENTER
    s = map(mapper.__getitem__, seq)
    s = ''.join(s)
    return s

def seq__reverseComplement(seq,mapper=None):
    seq = seq__complement(seq)
    seq = seq[::-1]
    return seq




#####################
# 20190812 splicing functions
# from HTSeq import GenomicInterval
# from pyRiboSeq.htseq_comp import iv__intersect, iv__null

_null = iv__null()
def iv__flip(iv):
    '''
    Note the GenomeInterval.__repr__() needs modification to reflect negative indexes
    [0,10) -> (-10,0], but repr [-10,0)
    by convention, we mentally map [-10,0) -> (10,0], but [-10,0) is still good for avoiding overlap
    '''
    assert iv.strand != '.',('Cannot flip a strandless interval',iv)
    iv = iv.copy()
    iv.start, iv.end = -iv.end, -iv.start
    iv.strand = {'+':'-','-':'+'}[iv.strand]    
    return iv

def iv__iv__split(ivq,ivr):
    '''
    Purpose:Split ivq(:GenomicInterval) to return (left_overhang, middle, right_overhang)
    
    Note (see also testing):
    no-overhang
    ivq:----xxxx-----
    ivr:----xxxxxx---
    
    left-overhang:
    ivq:--xxxxx------
    ivr:----xxxxxx---
    
    right-overhang:
    ivq:----xxxxxxx--
    ivr:----xxxxxx---
    
    left-right-overhang:
    ivq:--xxxxxxxxx--
    ivr:----xxxxxx---

    '''
    _null = iv__null()    
    if ivq.chrom != ivr.chrom or ivq.strand != ivr.strand:
        assert 0,(ivq,ivr)
    else:
#         if ivq.chrom=="-":
#             ivq,ivr = map(iv__flip,[ivq,ivr])
        chrom = ivq.chrom
        strand = ivq.strand
        ### get left
        if ivq.start < ivr.start:
            _left = GenomicInterval( chrom, ivq.start, ivr.start, strand)
        else:
            _left = _null
        
        ### get middle
        _start = max(ivr.start, ivq.start)
        _end   = min(ivr.end, ivq.end)
        if _start < _end:
            _mid = GenomicInterval(chrom, _start,_end, strand)
        else:
            _mid = _null
        
        ### get right
        if ivq.end > ivr.end:
            _right = GenomicInterval(chrom, ivr.end, ivq.end, strand)
        else:
            _right = _null
        return (_left,_mid,_right)

if __name__ == '__main__':
    # from pymisca.ext import ppJson
    import json
    def dict__assert(d):
        assert d['func']()==d['exp'], dppJson([d['func'](),d['exp']])
    def ppJson(d,**kw):
        '''
        Pretty print a dictionary
        '''
        s = json.dumps(d,indent=4, sort_keys=True,**kw)
        return s
    def dppJson(d,default=repr,**k):
        return ppJson(d,default=default,**k)    
    
    ####### testing: iv__iv__split()
    ### no overhang
    ia = GenomicInterval("NA",0,10)        
    ib = GenomicInterval("NA",0,10)
    exp = (_null, ia, _null,)
    res = iv__iv__split(ia,ib)
    assert res==exp,ppJson([res,exp],default=repr)

    ### no overhang
    ia = GenomicInterval("NA",2,8)        
    ib = GenomicInterval("NA",0,10)
    exp = (_null, ia, _null,)
    res = iv__iv__split(ia,ib)
    assert res==exp,ppJson([res,exp],default=repr)

    ### left overhang
    ia = GenomicInterval("NA",-1,10)        
    ib = GenomicInterval("NA",0,10)
    exp = ( GenomicInterval("NA",-1,0), 
           ib, _null,)
    res = iv__iv__split(ia,ib)
    assert res==exp,ppJson([res,exp],default=repr)


    ### right-overhang
    ia = GenomicInterval("NA",0,11)        
    ib = GenomicInterval("NA",0,10)
    exp = ( GenomicInterval("NA",0,0), 
           ib,  
           GenomicInterval("NA",10,11),
          )
    res = iv__iv__split(ia,ib)
    assert res==exp,ppJson([res,exp],default=repr)


    ### left-right-overhang
    ia = GenomicInterval("NA",-1,11)        
    ib = GenomicInterval("NA",0,10)
    exp = ( GenomicInterval("NA",-1,0), 
            ib,  
            GenomicInterval("NA",10,11),
          )
    res = iv__iv__split(ia,ib)
    assert res==exp,ppJson([res,exp],default=repr)

    #### test flipping
    ia = GenomicInterval("NA",-1,11,'-')
    ib = GenomicInterval("NA",0,10,'-')
    exp = ( GenomicInterval("NA",-1,0,'-'), 
            ib,  
            GenomicInterval("NA",10,11,'-'),
           )

    res = iv__iv__split(ia,ib)
    assert res==exp,ppJson([res,exp,],default=repr)
    
    #### aka  result becomes directional if wrapped in iv__flip
    
    res = tuple(map(iv__flip,iv__iv__split(*map(iv__flip,(ia,ib)))))
    assert res==exp[::-1],ppJson([res,exp[::-1],],default=repr)