#!/usr/bin/env bash
ml metal
metal ./type1.metal >../result/metal_type1.log &
metal ./type2.metal >../result/metal_type2.log &
metal ./type3.metal >../result/metal_type3.log &
metal ./type4.metal >../result/metal_type4.log &
metal ./type5.metal >../result/metal_type5.log &