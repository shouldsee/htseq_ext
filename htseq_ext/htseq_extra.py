#!/usr/bin/env python2

# pyext.execBaseFile('headers/header__import.py')

# if 1:
# stranded= 'yes'
import sys,os,re,glob,copy
import itertools,collections,functools
import HTSeq
# import pyRiboSeq.htseq_comp
import htseq_ext.htseq_comp as htseq_comp
# htseq_comp = pyRiboSeq.htseq_comp
# import copy
import funcy

import warnings
from HTSeq import GenomicInterval

class GenomicIntervalDeque(collections.deque):
    '''
    An ordered iterable of HTSeq.GenomicInterval(), intended to represent a transcript.
    
    property:
        cigarstring: inherited from an alignment if applicable
        iv: merge contained intervals to a big GenomicInterval()
        
    method:
        defaults inherited from collections.deque 
        fromCigar: init from parsing a cigar string 
        fromAlignment: init from parsing pysam.AlignedSegment().cigarstring
        countMissWithModel: assuming self is a read, count mismatch against a model
        countMissWithRead: assuming self is a model, ccount mismatch against a read
        (deprecated) compatibleWithModel: assuming self is a read, check compatibility against a model
        (deprecated) compatibleWithRead: assuming self is a model, check compatibility against a read
    
    
    '''
    def toTuples(self):
        res = []
        for x in self:
            v = (x.chrom,x.start,x.end,x.strand)
            res.append(v)
        return res
    @classmethod
    def fromTuples(cls,v):
        res = cls([])
        for x in v:
            vv = GenomicInterval(*x)
#             x.chrom,x.start,x.end,x.strand)
#             vv = GenomicInterval(x.chrom,x.start,x.end,x.strand)
            res.append(vv)
        return res

    ### +++ extra methods for Deque()
    def __getstate__(self):
        d = dict(self.__dict__)
        return d
    
    def __setstate__(self, d):
        ''' unpickle
         #### I *think* this is a safe way to do it
         ## This has not broken anything yet
        '''
        self.__dict__.update(d)
        
    def copy(self,):
        #### [FRAG] this is a shallow copy
        res = copy.copy(self)
        return res
    
    def sort(self, *args, **kwargs): 
        it = (self.pop() for x in xrange(len(self)))
        self.extend(sorted(it, *args, **kwargs))
    ### --- extra methods for Deque()
        
        
    def __init__(self, 
                 iterable,
                 maxlen=None,
                 name = None,
                 cigarstring = None,
                 sortKey= None,                 
                 strand = None,
                 last = None,
                 **kw):
        super(GenomicIntervalDeque,self,).__init__(iterable=iterable,maxlen=maxlen)
#         if not 
        if name is None:
            name = 'test'
        self.name = name
        self.cigarstring = cigarstring
        self.last = last
        
        self.strand = strand
        for iv in self:
            self._check(iv)
            
        if sortKey is not None:
            self.sortBy(sortKey)
            
    def _check(self, iv):
        if self.strand is None:
            self.strand = iv.strand
        else:
            assert iv.strand == self.strand,'interval {iv} should have strand:{strand}'.format(**locals())
         
    ####
    @property
    def length(self):
        return sum(x.length for x in self)
    @property
    def coordSum(self):
        if (self.end - self.start) != self.length:
            warnings.warn('coordSum evaulated for a non-contiguous GenomicIntervalDeque():%s'%self)
        return (self.start + self.end)
    @property
    def span(self):
        return (self.start,self.end)
    
    #### Start and end
    @property
    def startInterval(self):
        return self[0]    
    @property
    def startInterval_d(self):
        if self.strand=='-':
            res = self[-1]
        else:
            res = self[0]
        return res
    
    @property
    def start(self):
        return self.startInterval.start
    @property
    def start_d_as_pos(self):
        return self.startInterval_d.start_d_as_pos
    
    @property
    def endInterval(self):
        res = self[-1]
        return res
    @property
    def endInterval_d(self):
        if self.strand=='-':
            res = self[0]
        else:
            res = self[-1]
        return res
    @property
    def end(self):
        return self.endInterval.end
    @property
    def end_d_as_pos(self):
        return self.endInterval_d.end_d_as_pos


        
    def sortBy(self, key = None):
        if key is None:
            key = lambda x:x.start
        
        self.sort(key=key)
        return self
    
    def leftIter(self, 
                 strand='+',
                 **kwargs
                ):
        return self.getStrandedIter(strand=strand,**kwargs)
    
    def getStrandedIter(self, 
                 strand = None,
                 copy=True,
                 enum = True,
                 addSkipLen = True,
                ):
        if copy:
            _self =  self.copy()
        else:
            _self = self
        if strand is None:
            strand = self.strand
        assert strand is not None,'Must specify strand'
            
        def reducer(d,i):
#             if d is funcy.primitives.EMPTY:
#                 d = {}
            d = d.copy()
            if d == {}:
                d['iv'] = None
                d['skipLen'] = 0
            if enum:
                d['i'] = i
                
            if _self:
                if addSkipLen:
                    if d['iv'] is not None:
                        d['skipLen'] += d['iv'].length
#                     d['start'] += len()

                if htseq_comp.iv__isReverse(strand=strand):
                    d['iv'] = _self.pop()
                else:
                    ### assume positive strand
                    d['iv'] = _self.popleft()
                    

                ## [TODO] popright if negative strand
            else:
                d['iv'] = None
                
            return d
        
        it = funcy.ireductions(reducer, range(len(self)+1), acc = {} )
        
        return it
    
    @property
    def iv(self,):
        _iv = self[0].copy()
        _iv.extend_to_include(self[-1])
        return _iv
    @property
    def ivdq(self,):
        _self = self.copy()
        _self.clear()
        _self.append(self.iv)
        return _self
    
    @classmethod
    def fromCigar(cls, cigar, drop=None,dropEmpty=None, start=None,chrom=None,strand=None):
        if drop is None:
            drop = 'DNS'

        assert isinstance(cigar,basestring),cigar
        args = [cigar]
        if start is not None:
            args += [start]
            if chrom is not None:
                args += [chrom]
                if strand is not None:
                    args += [strand]
        cos = HTSeq.parse_cigar(*args)
        cos = [co for co in cos if co.type not in drop]
        
#         if dropEmpty is None:
#             dropEmpty = True
#         if dropEmpty:
#             cos = [co for co in cos if co.size ]  
        
        ivdq = GenomicIntervalDeque([co.ref_iv for co in cos ], cigarstring = cigar)
#         if start is not None:
#             for iv in ivdq:
#                 iv.start+= start
#                 iv.end  += start
        return ivdq
    
    @classmethod
    def fromAlignment(cls,aln, drop=None):
        ivdq = cls.fromCigar(aln.cigarstring, 
                                   start=aln.reference_start, 
                                   chrom=aln.reference_name, 
                                   strand=['+','-'][aln.is_reverse])    
        return ivdq
    def compatibleWithModel(self, ivdq_model,**kwargs):
        ivdq_read = self
        return htseq_comp.model__read__compatibleSpliced(ivdq_model, ivdq_read,**kwargs)
    
    def compatibleWithRead(self, ivdq_read, **kwargs):
        ivdq_model = self
        return htseq_comp.model__read__compatibleSpliced(ivdq_model, ivdq_read,**kwargs)
    
    def countMissWithModel(self,ivdq_model,**kwargs):
        ivdq_read = self
        res = htseq_comp.model__read__countMismatch(
            ivdq_model,ivdq_read,**kwargs
        )                                                  
        return res

    def countMissWithRead(self,ivdq_read,**kwargs):
        ivdq_model = self
        res = htseq_comp.model__read__countMismatch(
            ivdq_model,ivdq_read,**kwargs
        )                                                  
        return res
    
    def rebaseWithModel(self, ivdq_model, copy=False, **kwargs):
        '''
        Dont copy "self"
        '''
        ivdq_read = self
        res = htseq_comp.model__read__rebase( ivdq_model, ivdq_read, copy=copy, **kwargs)
        return res
    
    def rebaseWithRead(self, ivdq_read, **kwargs):
        ivdq_model = self
        res = htseq_comp.model__read__rebase( ivdq_model, ivdq_read, **kwargs)
        return res
        
    
###### 20190812
## rewrote splicing functions
## 

import collections
_DICT_CLASS = collections.OrderedDict
# from pyRiboSeq.htseq_extra import GenomicIntervalDeque
from htseq_ext.htseq_comp import iv__iv__split, iv__flip, _null

class ValuedIterator(object):
    _START = object()
    _END = object()
    def __init__(self,it):
        self.it = iter(it)
        self.value = self._START
        self.inserted = collections.deque()
    @property
    def ended(self):
        return self.value is self._END
    
    @property
    def started(self):
        return self.value is not self._START
    def next(self):
        return self.__next__()
    def insert(self,v):
        self.inserted.append(v)
        
    def __next__(self):
        if self.inserted:
            self.value = self.inserted.pop()
        else:
            self.value = value = next(self.it, self._END)
        return value
    
    

def ivdqq__ivdqr__rebase(ivqs,ivrs):
    return ivdqq__ivdqr__splicePairs(ivqs,ivrs)

def ivdqq__ivdqr__splice(ivqs,ivrs):
    res = ivdqq__ivdqr__splicePairs(ivqs,ivrs)
    res = ivPairs__splice(res)
    return GenomicIntervalDeque([x[0] for x in res])

def ivdqq__ivdqr__splicePairs(ivqs,ivrs):
    '''
    Filter out
    '''
    ivqit = ValuedIterator(ivqs)
    ivrit = ValuedIterator(ivrs)
    ivqit.next();
    ivrit.next();
#     out = GenomicIntervalDeque([])
    out = []
    # while True:
    while not ivrit.ended and not ivqit.ended:
        ivq = ivqit.value
        ivr = ivrit.value
        l,m,r = iv__iv__split(ivq,ivr)
        if m == _null:
            if r != _null:
                ivrit.next();
            elif l!= _null:
                ivqit.next();
            else:
                assert 0,dppJson((l,m,r,ivq,ivr))
        else:
            ### shift coordinate by ivr later
            out.append((m,ivr))
            if r != _null:
                ### avoiding ivqit.insert
                ivqit.value = r
                ivrit.next()
            else:
                ivqit.next()
    return out


def ivPairs__splice(ivPairs, inverse=False):
#     out = GenomicIntervalDeque([])
    out = []
    for i,(ivq,ivr) in enumerate(ivPairs):
        if not i:
            ivrl = ivr
            mapper = (0, ivr.start)
        elif ivr != ivrl:
            mapper = ( mapper[0] + ivrl.length, 
                      ivr.start)
#             mapper[0] += ivrl.length
#             mapper[1] = ivr.start 
        if not inverse:
            shifter = mapper[0] - mapper[1]
        else:
            shifter = mapper[1] - mapper[0]
            
        ivq = ivq.copy()
        ivq.start += shifter
        ivq.end += shifter
        out.append((ivq,ivr))
        ivrl = ivr        
    return out

def ivPairs__desplice(ivPairs):
    return ivPairs__splice(ivPairs,inverse=True)

def ivdqq__ivdqr__desplice(ivqs,ivqr):
    
    ### self splicing
    ivrs_dict = _DICT_CLASS(x for x in ivPairs__splice(zip(ivrs,ivrs)))
    ivrs_spliced = ivrs_dict.keys()
    
    ivout_pairs = ivdqq__ivdqr__splicePairs(ivqs, ivrs_spliced,)        
    ## mapping back to despliced "ivr"
    ivout_temp = [(ivq,ivrs_dict[ivr]) for ivq, ivr in ivout_pairs]
    
    ## desplice to recover ivqs
    ivout_recover = res = ivPairs__desplice(ivout_temp)
    return GenomicIntervalDeque([x[0] for x in ivout_recover])

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
        
    ### iv__rebase
    # ivrs = iva = GenomicIntervalDeque.fromCigar('20M5D10M',start=0,chrom="NA",strand="+")
    # ivqs = ivb = GenomicIntervalDeque.fromCigar('10M15D10M',start=0,chrom="NA",strand="+")
    # ivout = ivdqq__ivdqr__rebase(ivqs,ivrs)
    # print(dppJson([ivrs,ivqs,ivout]))        
    # ivrs = iva = GenomicIntervalDeque.fromCigar('20M5D10M',start=0,chrom="NA",strand="+")
    # ivqs = ivb = GenomicIntervalDeque.fromCigar('5D10M2D12M2D8M',start=0,chrom="NA",strand="+")
    ivrs = GenomicIntervalDeque.fromTuples(
        [('NA', 0, 20, '+'), 
         ('NA', 25, 35, '+')]
    )
    ivqs = GenomicIntervalDeque.fromTuples(
        [('NA', 5, 15, '+'), 
         ('NA', 17, 29, '+'), 
         ('NA', 31, 39, '+')]
    )

    ### expected coordinates after splicing against "ivrs"
    ### note there is loss of information
    ivqs_spliced = GenomicIntervalDeque.fromTuples(
    [('NA', 5, 15, '+'),
     ('NA', 17, 20, '+'),
     ('NA', 20, 24, '+'),
     ('NA', 26, 30, '+')]
    )

    ivout = ivdqq__ivdqr__splicePairs(ivqs,ivrs)
    # print(dppJson([ivrs,ivqs,ivout]))

    ivout_spliced  = res = ivPairs__splice(ivout)
    res = GenomicIntervalDeque([x[0] for x in res]).toTuples()
    exp = ivqs_spliced.toTuples()
    assert res == exp,dppJson((res,exp,))

    ivout_despliced  = res= ivPairs__desplice(ivout_spliced)
    exp = ivout
    assert res == exp,dppJson((res,exp,))


    "ivdqq__ivdqr__desplice"
    ivqs_recover = res = ivdqq__ivdqr__desplice(ivqs_spliced,ivrs)
    exp = GenomicIntervalDeque.fromTuples([('NA', 5, 15, '+'),
     ('NA', 17, 20, '+'),
     ('NA', 25, 29, '+'),
     ('NA', 31, 35, '+')])
    # assert 
    assert res == exp,dppJson((res,exp,))


    d = dict(func = lambda: ivdqq__ivdqr__splice(ivqs,ivrs),
         exp  = ivqs_spliced 
        )

    dict__assert(d)    
###
############

                
#### The context manager should process
class genomicSession(object):
    '''
    Manage the forward and backward transformation of the coordinates
    '''
    def __init__(self, ivdqs_model, ax=None, ):
        
        if isinstance(ivdqs_model, GenomicIntervalDeque):
            ### wrap single element
            ivdqs_model = [ivdqs_model]
        
        assert self.checkNames(ivdqs_model)
        self.transformers = self.ivdqs_model = ivdqs_model
    def __repr__(self):
        return unicode(self.transformers)

    
    
    @staticmethod
    def checkNames(ivdqs_model):
        for i, ivdq in enumerate(ivdqs_model):
            if i==0:
                res = True
            else:
                res &= (ivdq.chrom == ivdq_last.name)
            ivdq_last = ivdq
            if not res:
                break
                
        return res
        
    def __call__(self, ivdq,**kw):
        return self.forward(ivdq,**kw)
    
    def forward(self, ivdq,**kw):
        res = reduce(lambda ivdq,m: m.rebaseWithRead(ivdq, **kw),
                     self.ivdqs_model,
                     ivdq)
        return res
    
    def inverse(self,ivdq,**kw):
        raise Exception('Not implemented')
        
    def _filterKeep(self, ivdq, maxMiss, **kw):
        nMiss = 0
        _ivdq = ivdq
        if not _ivdq:
            return False
        
        for ivdq_model in self.transformers:
            nMiss += ivdq_model.countMissWithRead( _ivdq, maxMiss=maxMiss,**kw)
            if nMiss >= maxMiss:
                return False
            else:
                _ivdq = ivdq_model.rebaseWithRead(_ivdq)                
        return True
    
    def filterKeep(self, ivdqs, maxMiss = 5, **kw):
        ivdqs = copy.copy(ivdqs)
        
        res = filter(functools.partial( self._filterKeep, maxMiss=maxMiss, **kw),
                     ivdqs,
                    )
        return res

        
    def __enter__(self):
        return 
    def __exit__(self, type, value, traceback):
        return True  
    
    
    