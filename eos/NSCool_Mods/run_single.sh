#!/bin/sh

EOS=$1
TOV=$2
SF=$3

cp -r Model_template Model_1
sed -i -e "s/ENTER_EOS_HERE/$EOS/g" Model_1/Cool_Try.in
sed -i -e "s/ENTER_TOV_HERE/$TOV/g" Model_1/Cool_Try.in
sed -i -e "s/ENTER_SF_HERE/$SF/g" Model_1/Cool_Try.in

./NSCool.out

NAME=${EOS}_${TOV}_${SF}
rm -rf $NAME
mv Model_1 $NAME
