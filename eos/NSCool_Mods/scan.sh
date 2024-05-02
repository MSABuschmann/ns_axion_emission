#!/bin/sh

for SF in 0-0-0 SFB-0-0 SFB-0-T73 SFB-a-T73 SFB-b-T73 SFB-c-T73; do
	./run_single.sh APR_EOS_Cat APR_Cat_1.4 $SF
	./run_single.sh BSk22_test BSk22_1.4 $SF
	./run_single.sh BSk24_test BSk22_1.4 $SF
	./run_single.sh BSk25_test BSk22_1.4 $SF
	./run_single.sh BSk26_test BSk22_1.4 $SF
done
