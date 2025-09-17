## Testing Different Cutoff Methods for hintiso

**Pipeline Cmd**

```zsh
cd ../isoforms/gong     
   
python3 exp_cutoffmethods.py \
  /Users/gongchen/Code/isoforms/models/worm.splicemodel \
  /Users/gongchen/Code/isoforms/bin/hmm \
  --hmm_model /Users/gongchen/Code/isoforms/models \
  --test_dir /Users/gongchen/Code/isoforms/gong/test_data \
  --max_files 0 \
  --k 5 \
  --p 10 \
  --t -5.0 \
  --min_sites 5 \
  --min_intron 35 \
  --min_exon 25 \
  --flank 99 \
  --limit 100 \
  > ../gong/result/re_cutmethods.res 2> run.log
```

**Results**   

[Click me](../gong/result/re_cutmethods.res)

## Testing What's the proper Cutoff Methods for percentile method

[exp_perc.py](../gong/exp_perc.py)

**Pipelne Cmd**

``` zsh
cd ../isoform/gong

python3 exp_perc.py \
  /Users/gongchen/Code/isoforms/models/worm.splicemodel \
  /Users/gongchen/Code/isoforms/bin/hmm \
  --hmm_model /Users/gongchen/Code/isoforms/models \
  --test_dir /Users/gongchen/Code/isoforms/gong/test_data \
  --max_files 20 \
  --percentiles 10 20 30 \
  --min_intron 35 \
  --min_exon 25 \
  --flank 99 \
  --limit 100 \
  --score_percentile 75 \
  --summary ../gong/result/re_percentile.res.summary \
  > ../gong/result/re_percentile.res 2> run.log
```

**Results**

[Don't touch me](../gong/result/re_percentile.res)

## Find the Optimal percentile value for different length of gene

[exp_optperc.py](../gong/exp_optperc.py)

**Pipeline Cmd**  

```zsh
python3 exp_optperc.py \
  /Users/gongchen/Code/isoforms/models/worm.splicemodel \
  /Users/gongchen/Code/isoforms/bin/hmm \
  --hmm_model /Users/gongchen/Code/isoforms/models \
  --test_dir /Users/gongchen/Code/isoforms/gong/test_data/smallgenes \
  --max_files 0 \
  --target_quality 90 \
  --percentile_step 5 \
  --min_percentile 5 \
  --max_percentile 50 \
  --score_percentile 75 \
  --summary ../gong/result/re_optperc.res \
  > ../gong/data/optper.tsv 2> run.long
```

**Data**

**Results**

[I am not result](../gong/result/re_opper.res)
