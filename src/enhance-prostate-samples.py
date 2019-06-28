import OS
import sys
import cooler

# path to cool file
low_resolution_samples =

c = cooler.Cooler('/Users/neda/prostate-samples/PCa13266.multi-res.cool')
resolution = c.binsize
#mat = c.matrix(balance=True).fetch('chr5:10,000,000-15,000,000')
