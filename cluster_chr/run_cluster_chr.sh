#!/bin/bash

usage()
{
    echo "*------------------------------------------------------------------------------------------------*"
    echo "| The information of hic and contig were used to cluster chromosomes"
    echo "|"
    echo "| Usage: `basename $0` -o"
    echo "|      -g: subgraphs of ctgs"
    echo "|      -r: REs and length of ctgs"
    echo "|      -a: allele pair of ctgs base partig"
    echo "|      -l: hic links of ctgs"
    echo "|      -c: chromosome number of ctgs"
    echo "|      -n: output file prefix"
    echo "*------------------------------------------------------------------------------------------------*"
    exit 0
}
while getopts ':g:r:a:l:c:n:' OPT;do
    case $OPT in
        g)
            gfa="$OPTARG";;
        r)
            re="$OPTARG";;
        a)
            allele="$OPTARG";;
        l)
            link="$OPTARG";;
        c)
            chr="$OPTARG";;
        n)
            name="$OPTARG";;
        *)
            usage;;
    esac
done

if [ -z $gfa ] || [ -z $re ] || [ -z $allele ] || [ -z $link ] || [ -z $chr ] || [ -z $name ] ; then
	usage
fi

python /data/duwenjie/opt/anHiC/cluster_chr/allele_nei.py \
                        -r ${re} \
                        -a ${allele}

python /data/duwenjie/opt/anHiC/cluster_chr/pipeline.allele.py \
                        -g ${gfa} \
                        -l ${link} \
                        -r ${re} \
                        -a allele_louvain_nei.csv \
                        -n ${name}
                    
python /data/duwenjie/opt/anHiC/cluster_chr/trans_cluster.py \
                        -c ${name}.allele.cluster.expand.txt \
                        -s group_ctgs_N100.txt \
                        -n ${name}

python /data/duwenjie/opt/anHiC/cluster_chr/pipeline.chr.py \
                        -c ${name}.allele.cluster.expand.txt \
                        -l ${link}  \
                        -s group_ctgs_N100.txt \
                        -n ${name}

python /data/duwenjie/opt/multilevel_cluster.py  \
                        -c ${name}.allele.hic.csv \
                        -o ${name}.chr.cluster.txt \
                        -r 1.5 

min_r=0.1
max_r=5
expected_result=0
expected_result=${chr}

while :;do

    if echo "$max_r > $min_r" | bc;then

        mid_r=$(echo "scale=2; ($max_r + $min_r) / 2" | bc)

        echo "$mid_r"

        python /data/duwenjie/opt/multilevel_cluster.py  \
                            -c ${name}.allele.hic.csv \
                            -o ${name}.chr.cluster.txt \
                            -r ${mid_r}
        
        result=$(wc -l < "${name}.chr.cluster.txt" | awk '{print $1}')

        if [[ "${result}" -eq "${expected_result}" ]]; then
            echo "Find the parameters for louvain: $mid_r"
            break
        elif [[ $result -lt $expected_result ]]; then
            min_r=$mid_r
        else
            max_r=$mid_r
        fi

        (( iteration_count++ ))

        if (( iteration_count >= 20 )); then
            echo "Reached the maximum number of iterations, exit search..."
            break
        fi
    fi
done

if [ $(wc -l < "${name}.chr.cluster.txt") -eq $chr ]; then
    python /data/duwenjie/opt/anHiC/cluster_chr/trans_allele_cluster.py \
                            -c1 ${name}.chr.cluster.txt \
                            -c2 ${name}.allele.cluster.ctg.txt \
                            -n ${name} 
else
    echo "Cluster ERROR... louvain was unable to cluster to the specified number of clusters!"
fi

echo "done!"