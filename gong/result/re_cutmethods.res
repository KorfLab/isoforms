# Geniso cutoff methods performance test
# Parameters: min_intron=35 min_exon=25 flank=99 limit=100
# Method params: k=5 p=10% t=-5.0 min_sites=5
# Files processed: 50
#
# Average performance across all files:
# Method   Overlap%  Don_reduce%  Acc_reduce%  Iso_ratio%
# -------  --------  -----------  -----------  ----------
# gap          10.4        61.0         83.0         19.2
# smooth       12.8        66.3         77.6         22.5
# perc         76.8        35.7         26.8         94.9
# topk         17.2        69.4         68.4         30.0
# thresh      100.0         0.0          0.0        100.0
# adapt        17.2        69.4         68.4         30.0