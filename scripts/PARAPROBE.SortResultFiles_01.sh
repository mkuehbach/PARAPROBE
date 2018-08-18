#!/bin/bash

##Script to sort results of a PARAPROBE job according to SimID
##Markus Kuehbach, 2018/08/17, m.kuehbach at mpie.de

#user input
simid=$1

#automatic part
dirnm=$simid

#build directory structure
eval "mkdir -p $dirnm"
eval "mkdir -p $dirnm/DescrStats"
eval "mkdir -p $dirnm/Clustering"
eval "mkdir -p $dirnm/DescrStats/RDF"
eval "mkdir -p $dirnm/DescrStats/1NN"
eval "mkdir -p $dirnm/DescrStats/RIPK"
eval "mkdir -p $dirnm/DescrStats/KNN"
eval "mkdir -p $dirnm/DescrStats/MKNN"
#eval "mkdir -p $dirnm/Clustering/MaxSepSummary"
eval "mkdir -p $dirnm/Clustering/MaxSepSizeDistr"
eval "mkdir -p $dirnm/Additional"

#copy files over if exist
eval "mv *SimID.$simid.RDF.* $dirnm/DescrStats/RDF"
eval "mv *SimID.$simid.1NN.* $dirnm/DescrStats/1NN"
eval "mv *SimID.$simid.RIPK.* $dirnm/DescrStats/RIPK"
eval "mv *SimID.$simid.KNN.* $dirnm/DescrStats/KNN" ##really okay
eval "mv *SimID.$simid.MKNN* $dirnm/DescrStats/MKNN"
eval "mv *SimID.$simid.MaxSepDmax* $dirnm/Clustering"
eval "mv *SimID.$simid.MaxSepClustSzDistr.* $dirnm/Clustering/MaxSepSizeDistr"
eval "mv *SimID.$simid.*MyProfiling.csv $dirnm/Additional"
eval "mv *SimID.$simid.*KDTree*.csv $dirnm/Additional"
eval "mv *SimID.$simid.*SpatialDistribution*.csv $dirnm/Additional"
