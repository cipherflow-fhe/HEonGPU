// Copyright 2024-2025 Alişah Özcan
// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0
// Developer: Alişah Özcan

#include "util.cuh"
#include <utility>
#include <map>

/**
 * @company CipherFlow
 */
std::vector<Data64> smallPrimes = {
2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
	73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
	179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281,
	283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409,
	419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541,
	547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659,
	661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809,
	811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941,
	947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069,
	1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223,
	1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, 1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373,
	1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511,
	1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, 1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657,
	1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, 1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811,
	1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889, 1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987,
	1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053, 2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129,
	2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213, 2221, 2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287,
	2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357, 2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423,
	2437, 2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531, 2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617,
	2621, 2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687, 2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741,
	2749, 2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819, 2833, 2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903,
	2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999, 3001, 3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079,
	3083, 3089, 3109, 3119, 3121, 3137, 3163, 3167, 3169, 3181, 3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257,
	3259, 3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331, 3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413,
	3433, 3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511, 3517, 3527, 3529, 3533, 3539, 3541, 3547, 3557, 3559, 3571,
	3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643, 3659, 3671, 3673, 3677, 3691, 3697, 3701, 3709, 3719, 3727,
	3733, 3739, 3761, 3767, 3769, 3779, 3793, 3797, 3803, 3821, 3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907,
	3911, 3917, 3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989, 4001, 4003, 4007, 4013, 4019, 4021, 4027, 4049, 4051, 4057,
	4073, 4079, 4091, 4093, 4099, 4111, 4127, 4129, 4133, 4139, 4153, 4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229, 4231,
	4241, 4243, 4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297, 4327, 4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409,
	4421, 4423, 4441, 4447, 4451, 4457, 4463, 4481, 4483, 4493, 4507, 4513, 4517, 4519, 4523, 4547, 4549, 4561, 4567, 4583,
	4591, 4597, 4603, 4621, 4637, 4639, 4643, 4649, 4651, 4657, 4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733, 4751,
	4759, 4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831, 4861, 4871, 4877, 4889, 4903, 4909, 4919, 4931, 4933, 4937,
	4943, 4951, 4957, 4967, 4969, 4973, 4987, 4993, 4999, 5003, 5009, 5011, 5021, 5023, 5039, 5051, 5059, 5077, 5081, 5087,
	5099, 5101, 5107, 5113, 5119, 5147, 5153, 5167, 5171, 5179, 5189, 5197, 5209, 5227, 5231, 5233, 5237, 5261, 5273, 5279,
	5281, 5297, 5303, 5309, 5323, 5333, 5347, 5351, 5381, 5387, 5393, 5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441, 5443,
	5449, 5471, 5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521, 5527, 5531, 5557, 5563, 5569, 5573, 5581, 5591, 5623, 5639,
	5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683, 5689, 5693, 5701, 5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791,
	5801, 5807, 5813, 5821, 5827, 5839, 5843, 5849, 5851, 5857, 5861, 5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939,
	5953, 5981, 5987, 6007, 6011, 6029, 6037, 6043, 6047, 6053, 6067, 6073, 6079, 6089, 6091, 6101, 6113, 6121, 6131, 6133,
	6143, 6151, 6163, 6173, 6197, 6199, 6203, 6211, 6217, 6221, 6229, 6247, 6257, 6263, 6269, 6271, 6277, 6287, 6299, 6301,
	6311, 6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361, 6367, 6373, 6379, 6389, 6397, 6421, 6427, 6449, 6451, 6469, 6473,
	6481, 6491, 6521, 6529, 6547, 6551, 6553, 6563, 6569, 6571, 6577, 6581, 6599, 6607, 6619, 6637, 6653, 6659, 6661, 6673,
	6679, 6689, 6691, 6701, 6703, 6709, 6719, 6733, 6737, 6761, 6763, 6779, 6781, 6791, 6793, 6803, 6823, 6827, 6829, 6833,
	6841, 6857, 6863, 6869, 6871, 6883, 6899, 6907, 6911, 6917, 6947, 6949, 6959, 6961, 6967, 6971, 6977, 6983, 6991, 6997,
	7001, 7013, 7019, 7027, 7039, 7043, 7057, 7069, 7079, 7103, 7109, 7121, 7127, 7129, 7151, 7159, 7177, 7187, 7193, 7207,
	7211, 7213, 7219, 7229, 7237, 7243, 7247, 7253, 7283, 7297, 7307, 7309, 7321, 7331, 7333, 7349, 7351, 7369, 7393, 7411,
	7417, 7433, 7451, 7457, 7459, 7477, 7481, 7487, 7489, 7499, 7507, 7517, 7523, 7529, 7537, 7541, 7547, 7549, 7559, 7561,
	7573, 7577, 7583, 7589, 7591, 7603, 7607, 7621, 7639, 7643, 7649, 7669, 7673, 7681, 7687, 7691, 7699, 7703, 7717, 7723,
	7727, 7741, 7753, 7757, 7759, 7789, 7793, 7817, 7823, 7829, 7841, 7853, 7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919,
	7927, 7933, 7937, 7949, 7951, 7963, 7993, 8009, 8011, 8017, 8039, 8053, 8059, 8069, 8081, 8087, 8089, 8093, 8101, 8111,
	8117, 8123, 8147, 8161, 8167, 8171, 8179, 8191, 8209, 8219, 8221, 8231, 8233, 8237, 8243, 8263, 8269, 8273, 8287, 8291,
	8293, 8297, 8311, 8317, 8329, 8353, 8363, 8369, 8377, 8387, 8389, 8419, 8423, 8429, 8431, 8443, 8447, 8461, 8467, 8501,
	8513, 8521, 8527, 8537, 8539, 8543, 8563, 8573, 8581, 8597, 8599, 8609, 8623, 8627, 8629, 8641, 8647, 8663, 8669, 8677,
	8681, 8689, 8693, 8699, 8707, 8713, 8719, 8731, 8737, 8741, 8747, 8753, 8761, 8779, 8783, 8803, 8807, 8819, 8821, 8831,
	8837, 8839, 8849, 8861, 8863, 8867, 8887, 8893, 8923, 8929, 8933, 8941, 8951, 8963, 8969, 8971, 8999, 9001, 9007, 9011,
	9013, 9029, 9041, 9043, 9049, 9059, 9067, 9091, 9103, 9109, 9127, 9133, 9137, 9151, 9157, 9161, 9173, 9181, 9187, 9199,
	9203, 9209, 9221, 9227, 9239, 9241, 9257, 9277, 9281, 9283, 9293, 9311, 9319, 9323, 9337, 9341, 9343, 9349, 9371, 9377,
	9391, 9397, 9403, 9413, 9419, 9421, 9431, 9433, 9437, 9439, 9461, 9463, 9467, 9473, 9479, 9491, 9497, 9511, 9521, 9533,
	9539, 9547, 9551, 9587, 9601, 9613, 9619, 9623, 9629, 9631, 9643, 9649, 9661, 9677, 9679, 9689, 9697, 9719, 9721, 9733,
	9739, 9743, 9749, 9767, 9769, 9781, 9787, 9791, 9803, 9811, 9817, 9829, 9833, 9839, 9851, 9857, 9859, 9871, 9883, 9887,
	9901, 9907, 9923, 9929, 9931, 9941, 9949, 9967, 9973, 10007, 10009, 10037, 10039, 10061, 10067, 10069, 10079, 10091, 10093, 10099,
	10103, 10111, 10133, 10139, 10141, 10151, 10159, 10163, 10169, 10177, 10181, 10193, 10211, 10223, 10243, 10247, 10253, 10259, 10267, 10271,
	10273, 10289, 10301, 10303, 10313, 10321, 10331, 10333, 10337, 10343, 10357, 10369, 10391, 10399, 10427, 10429, 10433, 10453, 10457, 10459,
	10463, 10477, 10487, 10499, 10501, 10513, 10529, 10531, 10559, 10567, 10589, 10597, 10601, 10607, 10613, 10627, 10631, 10639, 10651, 10657,
	10663, 10667, 10687, 10691, 10709, 10711, 10723, 10729, 10733, 10739, 10753, 10771, 10781, 10789, 10799, 10831, 10837, 10847, 10853, 10859,
	10861, 10867, 10883, 10889, 10891, 10903, 10909, 10937, 10939, 10949, 10957, 10973, 10979, 10987, 10993, 11003, 11027, 11047, 11057, 11059,
	11069, 11071, 11083, 11087, 11093, 11113, 11117, 11119, 11131, 11149, 11159, 11161, 11171, 11173, 11177, 11197, 11213, 11239, 11243, 11251,
	11257, 11261, 11273, 11279, 11287, 11299, 11311, 11317, 11321, 11329, 11351, 11353, 11369, 11383, 11393, 11399, 11411, 11423, 11437, 11443,
	11447, 11467, 11471, 11483, 11489, 11491, 11497, 11503, 11519, 11527, 11549, 11551, 11579, 11587, 11593, 11597, 11617, 11621, 11633, 11657,
	11677, 11681, 11689, 11699, 11701, 11717, 11719, 11731, 11743, 11777, 11779, 11783, 11789, 11801, 11807, 11813, 11821, 11827, 11831, 11833,
	11839, 11863, 11867, 11887, 11897, 11903, 11909, 11923, 11927, 11933, 11939, 11941, 11953, 11959, 11969, 11971, 11981, 11987, 12007, 12011,
	12037, 12041, 12043, 12049, 12071, 12073, 12097, 12101, 12107, 12109, 12113, 12119, 12143, 12149, 12157, 12161, 12163, 12197, 12203, 12211,
	12227, 12239, 12241, 12251, 12253, 12263, 12269, 12277, 12281, 12289, 12301, 12323, 12329, 12343, 12347, 12373, 12377, 12379, 12391, 12401,
	12409, 12413, 12421, 12433, 12437, 12451, 12457, 12473, 12479, 12487, 12491, 12497, 12503, 12511, 12517, 12527, 12539, 12541, 12547, 12553,
	12569, 12577, 12583, 12589, 12601, 12611, 12613, 12619, 12637, 12641, 12647, 12653, 12659, 12671, 12689, 12697, 12703, 12713, 12721, 12739,
	12743, 12757, 12763, 12781, 12791, 12799, 12809, 12821, 12823, 12829, 12841, 12853, 12889, 12893, 12899, 12907, 12911, 12917, 12919, 12923,
	12941, 12953, 12959, 12967, 12973, 12979, 12983, 13001, 13003, 13007, 13009, 13033, 13037, 13043, 13049, 13063, 13093, 13099, 13103, 13109,
	13121, 13127, 13147, 13151, 13159, 13163, 13171, 13177, 13183, 13187, 13217, 13219, 13229, 13241, 13249, 13259, 13267, 13291, 13297, 13309,
	13313, 13327, 13331, 13337, 13339, 13367, 13381, 13397, 13399, 13411, 13417, 13421, 13441, 13451, 13457, 13463, 13469, 13477, 13487, 13499,
	13513, 13523, 13537, 13553, 13567, 13577, 13591, 13597, 13613, 13619, 13627, 13633, 13649, 13669, 13679, 13681, 13687, 13691, 13693, 13697,
	13709, 13711, 13721, 13723, 13729, 13751, 13757, 13759, 13763, 13781, 13789, 13799, 13807, 13829, 13831, 13841, 13859, 13873, 13877, 13879,
	13883, 13901, 13903, 13907, 13913, 13921, 13931, 13933, 13963, 13967, 13997, 13999, 14009, 14011, 14029, 14033, 14051, 14057, 14071, 14081,
	14083, 14087, 14107, 14143, 14149, 14153, 14159, 14173, 14177, 14197, 14207, 14221, 14243, 14249, 14251, 14281, 14293, 14303, 14321, 14323,
	14327, 14341, 14347, 14369, 14387, 14389, 14401, 14407, 14411, 14419, 14423, 14431, 14437, 14447, 14449, 14461, 14479, 14489, 14503, 14519,
	14533, 14537, 14543, 14549, 14551, 14557, 14561, 14563, 14591, 14593, 14621, 14627, 14629, 14633, 14639, 14653, 14657, 14669, 14683, 14699,
	14713, 14717, 14723, 14731, 14737, 14741, 14747, 14753, 14759, 14767, 14771, 14779, 14783, 14797, 14813, 14821, 14827, 14831, 14843, 14851,
	14867, 14869, 14879, 14887, 14891, 14897, 14923, 14929, 14939, 14947, 14951, 14957, 14969, 14983, 15013, 15017, 15031, 15053, 15061, 15073,
	15077, 15083, 15091, 15101, 15107, 15121, 15131, 15137, 15139, 15149, 15161, 15173, 15187, 15193, 15199, 15217, 15227, 15233, 15241, 15259,
	15263, 15269, 15271, 15277, 15287, 15289, 15299, 15307, 15313, 15319, 15329, 15331, 15349, 15359, 15361, 15373, 15377, 15383, 15391, 15401,
	15413, 15427, 15439, 15443, 15451, 15461, 15467, 15473, 15493, 15497, 15511, 15527, 15541, 15551, 15559, 15569, 15581, 15583, 15601, 15607,
	15619, 15629, 15641, 15643, 15647, 15649, 15661, 15667, 15671, 15679, 15683, 15727, 15731, 15733, 15737, 15739, 15749, 15761, 15767, 15773,
	15787, 15791, 15797, 15803, 15809, 15817, 15823, 15859, 15877, 15881, 15887, 15889, 15901, 15907, 15913, 15919, 15923, 15937, 15959, 15971,
	15973, 15991, 16001, 16007, 16033, 16057, 16061, 16063, 16067, 16069, 16073, 16087, 16091, 16097, 16103, 16111, 16127, 16139, 16141, 16183,
	16187, 16189, 16193, 16217, 16223, 16229, 16231, 16249, 16253, 16267, 16273, 16301, 16319, 16333, 16339, 16349, 16361, 16363, 16369, 16381,
	16411, 16417, 16421, 16427, 16433, 16447, 16451, 16453, 16477, 16481, 16487, 16493, 16519, 16529, 16547, 16553, 16561, 16567, 16573, 16603,
	16607, 16619, 16631, 16633, 16649, 16651, 16657, 16661, 16673, 16691, 16693, 16699, 16703, 16729, 16741, 16747, 16759, 16763, 16787, 16811,
	16823, 16829, 16831, 16843, 16871, 16879, 16883, 16889, 16901, 16903, 16921, 16927, 16931, 16937, 16943, 16963, 16979, 16981, 16987, 16993,
	17011, 17021, 17027, 17029, 17033, 17041, 17047, 17053, 17077, 17093, 17099, 17107, 17117, 17123, 17137, 17159, 17167, 17183, 17189, 17191,
	17203, 17207, 17209, 17231, 17239, 17257, 17291, 17293, 17299, 17317, 17321, 17327, 17333, 17341, 17351, 17359, 17377, 17383, 17387, 17389,
};

namespace heongpu
{
    bool coefficient_validator(const std::vector<int>& log_Q_bases_bit_sizes,
                               const std::vector<int>& log_P_bases_bit_sizes)
    {
        int P_size = log_P_bases_bit_sizes.size();

        int total_P_modulus = 0; // TODO: calculate it with prod
        for (int i = 0; i < P_size; i++)
        {
            total_P_modulus = total_P_modulus + log_P_bases_bit_sizes[i];
        }

        int Q_size = log_Q_bases_bit_sizes.size();

        int remainder = Q_size % P_size;
        int quotient = Q_size / P_size;

        int counter = 0;
        for (int i = 0; i < quotient; i++)
        {
            int pair_sum = 0;
            for (int j = 0; j < P_size; j++)
            {
                pair_sum = pair_sum + log_Q_bases_bit_sizes[counter];
                counter++;
            }

            /**
             * @company CipherFlow
             */
            if (pair_sum > 2 + total_P_modulus)
            {
                return false;
            }
        }

        int pair_sum = 0;
        for (int j = 0; j < remainder; j++)
        {
            pair_sum = pair_sum + log_Q_bases_bit_sizes[counter];
            counter++;
        }

        if (pair_sum > total_P_modulus)
        {
            return false;
        }

        return true;
    }

    Data64 extendedGCD(Data64 a, Data64 b, Data64& x, Data64& y)
    {
        if (a == 0)
        {
            x = 0;
            y = 1;
            return b;
        }

        Data64 x1, y1;
        Data64 gcd = extendedGCD(b % a, a, x1, y1);

        x = y1 - (b / a) * x1;
        y = x1;

        return gcd;
    }

    Data64 modInverse(Data64 a, Data64 m)
    {
        Data64 x, y;
        Data64 gcd = extendedGCD(a, m, x, y);

        if (gcd != 1)
        {
            // Modular inverse does not exist
            return 0;
        }
        else
        {
            // Ensure the result is positive
            Data64 result = (x % m + m) % m;
            return result;
        }
    }

    int calculate_bit_size(Data64 input)
    {
        return (int) log2(input) + 1;
    }

    bool is_power_of_two(size_t number)
    {
        return (number > 0) && ((number & (number - 1)) == 0);
    }

    int calculate_bit_count(Data64 number)
    {
        return log2(number) + 1;
    }

    int calculate_big_integer_bit_count(Data64* number, int word_count)
    {
        int size = word_count;
        for (int i = (word_count - 1); i > (-1); i--)
        {
            if (number[i] == 0)
            {
                size--;
            }
            else
            {
                break;
            }
        }

        return ((size - 1) * 64) + calculate_bit_count(number[size - 1]);
    }

    bool miller_rabin(const Data64& value, size_t num_rounds)
    {
        Modulus64 modulus_in(value);

        Data64 d = value - 1;
        Data64 r = 0;
        while (0 == (d & 0x1)) // #true while the last bit of r is zero
        {
            d >>= 1;
            r++;
        }
        if (r == 0)
        {
            return false;
        }

        // apply miller_rabin primality test
        std::random_device rand;
        std::uniform_int_distribution<Data64> dist(3, value - 1);
        for (size_t i = 0; i < num_rounds; i++)
        {
            Data64 a = i ? dist(rand) : 2;
            Data64 x = OPERATOR64::exp(a, d, modulus_in);
            if (x == 1 || x == value - 1)
            {
                continue;
            }
            Data64 count = 0;
            do
            {
                x = OPERATOR64::mult(x, x, modulus_in);
                count++;
            } while (x != value - 1 && count < r - 1);
            if (x != value - 1)
            {
                return false;
            }
        }
        return true;
    }

    bool is_prime(const Data64& value)
    {
        size_t num_rounds = 11;

        // First check the prime under 1000.
        std::vector<Data64> low_primes = {
            3ULL,   5ULL,   7ULL,   11ULL,  13ULL,  17ULL,  19ULL,  23ULL,
            29ULL,  31ULL,  37ULL,  41ULL,  43ULL,  47ULL,  53ULL,  59ULL,
            61ULL,  67ULL,  71ULL,  73ULL,  79ULL,  83ULL,  89ULL,  97ULL,
            101ULL, 103ULL, 107ULL, 109ULL, 113ULL, 127ULL, 131ULL, 137ULL,
            139ULL, 149ULL, 151ULL, 157ULL, 163ULL, 167ULL, 173ULL, 179ULL,
            181ULL, 191ULL, 193ULL, 197ULL, 199ULL, 211ULL, 223ULL, 227ULL,
            229ULL, 233ULL, 239ULL, 241ULL, 251ULL, 257ULL, 263ULL, 269ULL,
            271ULL, 277ULL, 281ULL, 283ULL, 293ULL, 307ULL, 311ULL, 313ULL,
            317ULL, 331ULL, 337ULL, 347ULL, 349ULL, 353ULL, 359ULL, 367ULL,
            373ULL, 379ULL, 383ULL, 389ULL, 397ULL, 401ULL, 409ULL, 419ULL,
            421ULL, 431ULL, 433ULL, 439ULL, 443ULL, 449ULL, 457ULL, 461ULL,
            463ULL, 467ULL, 479ULL, 487ULL, 491ULL, 499ULL, 503ULL, 509ULL,
            521ULL, 523ULL, 541ULL, 547ULL, 557ULL, 563ULL, 569ULL, 571ULL,
            577ULL, 587ULL, 593ULL, 599ULL, 601ULL, 607ULL, 613ULL, 617ULL,
            619ULL, 631ULL, 641ULL, 643ULL, 647ULL, 653ULL, 659ULL, 661ULL,
            673ULL, 677ULL, 683ULL, 691ULL, 701ULL, 709ULL, 719ULL, 727ULL,
            733ULL, 739ULL, 743ULL, 751ULL, 757ULL, 761ULL, 769ULL, 773ULL,
            787ULL, 797ULL, 809ULL, 811ULL, 821ULL, 823ULL, 827ULL, 829ULL,
            839ULL, 853ULL, 857ULL, 859ULL, 863ULL, 877ULL, 881ULL, 883ULL,
            887ULL, 907ULL, 911ULL, 919ULL, 929ULL, 937ULL, 941ULL, 947ULL,
            953ULL, 967ULL, 971ULL, 977ULL, 983ULL, 991ULL, 997ULL};

        if (value >= 3)
        {
            if ((value & 0x1) != 0)
            {
                for (int i = 0; i < low_primes.size(); i++)
                {
                    if (value == low_primes[i])
                    {
                        return true;
                    }
                    if ((value % low_primes[i]) == 0)
                    {
                        return false;
                    }
                }

                return miller_rabin(value, num_rounds);
            }
        }

        return false;
    }

    std::vector<Data64> generate_proper_primes(Data64 factor, int bit_size,
                                               size_t count)
    {
        std::vector<Data64> destination;

        // Start with (2^bit_size - 1) / factor * factor + 1
        Data64 value = ((Data64(0x1) << bit_size) - 1) / factor * factor + 1;

        Data64 lower_bound = Data64(0x1) << (bit_size - 1);
        while (count > 0 && value > lower_bound)
        {
            if (is_prime(value))
            {
                destination.emplace_back(std::move(value));
                count--;
            }
            value -= factor;
        }
        if (count > 0)
        {
            throw std::logic_error("failed to find enough qualifying primes");
        }
        return destination;
    }

    std::vector<Modulus64>
    generate_primes(size_t poly_modulus_degree,
                    const std::vector<int> prime_bit_sizes)
    {
        std::vector<Modulus64> prime_vector_;
        std::unordered_map<int, size_t> count_table;
        std::unordered_map<int, std::vector<Data64>> prime_table;
        for (int size : prime_bit_sizes)
        {
            if ((size > MAX_USER_DEFINED_MOD_BIT_COUNT) ||
                (size < MIN_USER_DEFINED_MOD_BIT_COUNT))
            {
                throw std::logic_error("invalid modulus bit size");
            }

            ++count_table[size];
        }

        Data64 factor = Data64(2) * Data64(poly_modulus_degree);
        for (const auto& table_elt : count_table)
        {
            prime_table[table_elt.first] = generate_proper_primes(
                factor, table_elt.first, table_elt.second);
        }

        for (int size : prime_bit_sizes)
        {
            prime_vector_.emplace_back(Modulus64(prime_table[size].back()));
            prime_table[size].pop_back();
        }

        return prime_vector_;
    }

    std::vector<Modulus64> generate_internal_primes(size_t poly_modulus_degree,
                                                    const int prime_count)
    {
        std::vector<Modulus64> all_primes;

        std::vector<int> prime_bit_sizes;
        for (int i = 0; i < prime_count; i++)
        {
            prime_bit_sizes.push_back(MAX_MOD_BIT_COUNT);
        }

        std::unordered_map<int, size_t> count_table;
        std::unordered_map<int, std::vector<Data64>> prime_table;
        for (int size : prime_bit_sizes)
        {
            ++count_table[size];
        }

        Data64 factor = Data64(2) * Data64(poly_modulus_degree);
        for (const auto& table_elt : count_table)
        {
            prime_table[table_elt.first] = generate_proper_primes(
                factor, table_elt.first, table_elt.second);
        }

        for (int size : prime_bit_sizes)
        {
            all_primes.emplace_back(Modulus64(prime_table[size].back()));
            prime_table[size].pop_back();
        }

        return all_primes;
    }

    /**
     * @company CipherFlow
     */
    Data64 gcd(Data64 a, Data64 b) {
        if (a == 0 || b == 0) {
            return 0;
        }

        while (b != 0)
        {
            a = a % b;
            std::swap(a, b);
        }

        return a;
        
    }

    /**
     * @company CipherFlow
     */
    Data64 polynomialPollardsRho(Data64 x1, Data64 x2, Data64 c) 
    {
        Data64 z = OPERATOR64::exp(x1, 2, x2);
        z += c;
        z %= x2;
        return z;
    }

    /**
     * @company CipherFlow
     */
    Data64 factorizationPollardsRho(Data64 m)
    {
        Data64 x, y, c, d;
        for (c = 1; c < 10; c++) {
            x = 2;
            y = 2;
            d = 1;

            while (d != 0)
            {
                x = polynomialPollardsRho(x, m, c);
                y = polynomialPollardsRho(polynomialPollardsRho(y, m, c), m, c);

                if (y > x) {
                    std::swap(x, y);
                }

                d = gcd(x-y, m);

                if (d > 1) {
                    return d;
                }
            }
            
        }

        return d;
    }

    /**
     * @company CipherFlow
     */
    std::vector<Data64> getFactors(Data64 n) {
        std::vector<Data64> factors;
        Data64 factor, m;
        m = n;

        for (int i = 0; i < smallPrimes.size(); i++) {
            Data64 smallPrime = smallPrimes[i];
            bool addFactor = false;
            while (m%smallPrime == 0)
            {
                m /= smallPrime;
			    addFactor = true;
            }

            if (addFactor) {
                factors.push_back(smallPrime);
            }
            
        }

        if (m==1) {
            return factors;
        }

        while (true)
        {
            factor = factorizationPollardsRho(m);
            if (factor == 0) {
                factors.push_back(m);
                break;
            }
            m /= factor;
            if (factors.size() > 0 && factor == factors.back()) {
                continue;
            }
            factors.push_back(factor);
        }
        
        return factors;
    }

    /**
     * @company CipherFlow
     */
    Data64 primitiveRoot(Data64 q) {
        Data64 tmp, g;

        bool notFoundPrimitiveRoot = true;

        std::vector<Data64> factors = getFactors(q - 1);

        g = 2;

        while (notFoundPrimitiveRoot)
        {
            g++;
            for (auto& factor : factors) {
                tmp = (q - 1) / factor;
                if (OPERATOR64::exp(g, tmp, q)==1) {
                    notFoundPrimitiveRoot = true;
				    break;
                }
                notFoundPrimitiveRoot = false;
            }
        }
        
        return g;
    }


    bool is_primitive_root(Data64 root, size_t degree, Modulus64& modulus)
    {
        // root^(degree/2) = modulus - 1 .
        Data64 degree_over2 = degree >> 1;

        return OPERATOR64::exp(root, degree_over2, modulus) ==
               (modulus.value - 1);
    }

    bool find_primitive_root(size_t degree, Modulus64& modulus,
                             Data64& destination)
    {
        Data64 size_entire_group = modulus.value - 1;

        Data64 size_quotient_group = size_entire_group / degree;

        if (size_entire_group - size_quotient_group * degree != 0)
        {
            return false;
        }

        std::random_device rd;

        int attempt_counter = 0;
        int attempt_counter_max = 100;
        do
        {
            attempt_counter++;

            Data64 random_num =
                (static_cast<Data64>(rd()) << 32) | static_cast<Data64>(rd());
            // destination = OPERATOR64::reduce(random_num, modulus);
            destination = random_num % modulus.value;

            // Raise the random number to power the size of the quotient
            // to get rid of irrelevant part
            destination =
                OPERATOR64::exp(destination, size_quotient_group, modulus);
        } while (!is_primitive_root(destination, degree, modulus) &&
                 (attempt_counter < attempt_counter_max));

        return is_primitive_root(destination, degree, modulus);
    }

    Data64 find_minimal_primitive_root(size_t degree, Modulus64& modulus)
    {
        Data64 root;
        if (!find_primitive_root(degree, modulus, root))
        {
            throw std::logic_error("no sufficient root unity");
        }

        Data64 generator_sq = OPERATOR64::mult(root, root, modulus);

        Data64 current_generator = root;

        for (size_t i = 0; i < degree; i += 2)
        {
            if (current_generator < root)
            {
                root = current_generator;
            }

            current_generator =
                OPERATOR64::mult(current_generator, generator_sq, modulus);
        }

        return root;
    }

    std::vector<Data64>
    generate_primitive_root_of_unity(size_t poly_modulus_degree,
                                     std::vector<Modulus64> primes)
    {
        std::vector<Data64> root_of_unity;

        // 2nth root of unity
        for (int i = 0; i < primes.size(); i++)
        {
            // root_of_unity.push_back(find_minimal_primitive_root(
            //     2 * poly_modulus_degree, primes[i]));
            
            /**
             * @company CipherFlow
             */
            Data64 g = primitiveRoot(primes[i].value);
            Data64 power = (primes[i].value - 1) / (2 * poly_modulus_degree);

            root_of_unity.push_back(OPERATOR64::exp(g, power, primes[i]));
        }

        return root_of_unity;
    }

    std::vector<Root64> generate_ntt_table(std::vector<Data64> psi,
                                           std::vector<Modulus64> primes,
                                           int n_power)
    {
        int n_ = 1 << n_power;
        std::vector<Root64> forward_table; // bit reverse order
        for (int i = 0; i < primes.size(); i++)
        {
            std::vector<Root64> table;
            table.push_back(1);

            for (int j = 1; j < n_; j++)
            {
                Data64 exp =
                    OPERATOR64::mult(table[(j - 1)], psi[i], primes[i]);
                table.push_back(exp);
            }

            for (int j = 0; j < n_; j++) // take bit reverse order
            {
                forward_table.push_back(table[gpuntt::bitreverse(j, n_power)]);
            }
        }

        return forward_table;
    }

    std::vector<Root64> generate_intt_table(std::vector<Data64> psi,
                                            std::vector<Modulus64> primes,
                                            int n_power)
    {
        int n_ = 1 << n_power;
        std::vector<Root64> inverse_table; // bit reverse order
        for (int i = 0; i < primes.size(); i++)
        {
            std::vector<Root64> table;
            table.push_back(1);

            Data64 inv_root = OPERATOR64::modinv(psi[i], primes[i]);
            for (int j = 1; j < n_; j++)
            {
                Data64 exp =
                    OPERATOR64::mult(table[(j - 1)], inv_root, primes[i]);
                table.push_back(exp);
            }

            for (int j = 0; j < n_; j++) // take bit reverse order
            {
                inverse_table.push_back(table[gpuntt::bitreverse(j, n_power)]);
            }
        }

        return inverse_table;
    }

    std::vector<Ninverse64> generate_n_inverse(size_t poly_modulus_degree,
                                               std::vector<Modulus64> primes)
    {
        Data64 n_ = poly_modulus_degree;
        std::vector<Ninverse64> n_inverse;
        for (int i = 0; i < primes.size(); i++)
        {
            n_inverse.push_back(OPERATOR64::modinv(n_, primes[i]));
        }

        return n_inverse;
    }

    __global__ void unsigned_signed_convert(Data64* input, Data64* output,
                                            Modulus64* modulus)
    {
        int idx = blockIdx.x * blockDim.x + threadIdx.x; // ring size

        int64_t threshold = (modulus[0].value + 1) >> 1;
        int64_t input_reg = static_cast<int64_t>(input[idx]);
        input_reg = (input_reg > threshold)
                        ? input_reg - static_cast<int64_t>(modulus[0].value)
                        : input_reg;

        output[idx] = static_cast<Data64>(input_reg);
    }

    __global__ void fill_device_vector(Data64* vector, Data64 number, int size)
    {
        int idx = blockIdx.x * blockDim.x + threadIdx.x; // ring size

        if (idx < size)
        {
            vector[idx] = number;
        }
    }

    int find_closest_divisor(int N)
    {
        double target = std::sqrt(N);
        int closest_div = 1;
        double min_diff = std::abs(closest_div - target);

        for (int k = 1; k <= std::sqrt(N); ++k)
        {
            if (N % k == 0)
            {
                for (int divisor : {k, N / k})
                {
                    double diff = std::abs(divisor - target);
                    if (diff < min_diff)
                    {
                        min_diff = diff;
                        closest_div = divisor;
                    }
                }
            }
        }
        return closest_div;
    }

    std::vector<std::vector<int>> split_array(const std::vector<int>& array,
                                              int chunk_size)
    {
        std::vector<std::vector<int>> result;
        int n = array.size();
        for (int i = 0; i < n; i += chunk_size)
        {
            result.emplace_back(array.begin() + i,
                                array.begin() + min(i + chunk_size, n));
        }
        return result;
    }

    std::vector<std::vector<int>> seperate_func(const std::vector<int>& A)
    {
        int initial_size = A.size();
        int counter = 2;
        int offset = A[1] - A[0];

        for (size_t i = 1; i < A.size() - 1; ++i)
        {
            if (A[i + 1] - A[i] != offset)
            {
                break;
            }
            counter++;
        }

        int real_n1 = heongpu::find_closest_divisor(counter);

        if (counter == initial_size)
        {
            return split_array(A, real_n1);
        }
        else
        {
            auto first_part = split_array(
                std::vector<int>(A.begin(), A.begin() + counter), real_n1);
            auto second_part = split_array(
                std::vector<int>(A.begin() + counter, A.end()), real_n1);

            first_part.insert(first_part.end(), second_part.begin(),
                              second_part.end());
            return first_part;
        }
    }

    std::vector<std::vector<int>> bsgs_index(const std::vector<int>& array,
                                             int N1, std::vector<int>& rot_n1,
                                             std::vector<int>& rot_n2)
    {
        std::vector<std::vector<int>> result;
        int n = array.size();
       
        std::map<int, bool> rot_n1_map;
        std::map<int, bool> rot_n2_map;

        std::vector<int> temp;
        int prev_idx_n1 = (array[0] / N1) * N1;
        for (const auto& rot : array)
        {   
            int idx_n1 = (rot / N1) * N1;
            if (idx_n1 != prev_idx_n1)
            {
                result.push_back(temp);
                temp.clear();
                prev_idx_n1 = idx_n1;
            }

            int idx_n2 = rot % N1;

            temp.push_back(idx_n1 + idx_n2);
           
            rot_n1_map[idx_n1] = true;
            rot_n2_map[idx_n2] = true;
        }

        result.push_back(temp);

        for (const auto& [key, _] : rot_n1_map)
        {
            rot_n1.push_back(key);
        }
        for (const auto& [key, _] : rot_n2_map)
        {
            rot_n2.push_back(key);
        }
        return result;
    }

int find_best_bsgs_split(const std::vector<int>& array, int max_N,
                             float bsgs_ratio)
    {
        const float epsilon = 1e-8f;  

        for (int N1 = 1; N1 < max_N; N1 <<= 1)
        {
            std::vector<int> rot_n1;
            std::vector<int> rot_n2;
            bsgs_index(array, N1, rot_n1, rot_n2);

            float current_ratio = float(rot_n2.size() - 1) / float(rot_n1.size() - 1);

            if (std::abs(current_ratio - bsgs_ratio) < epsilon)
            {
                return N1;
            }
            if (current_ratio > bsgs_ratio)
            {
                return N1 / 2;
            }
        }

        return 1;
    }

    std::vector<std::vector<int>> seperate_func_v2(const std::vector<int>& A,
                                                   int slots,
                                                   std::vector<int>& rot_n1,
                                                   std::vector<int>& rot_n2,
                                                   float bsgs_ratio)
    {
        int N1 = find_best_bsgs_split(A, slots, bsgs_ratio);
        return bsgs_index(A, N1, rot_n1, rot_n2);
    }

    std::vector<int> unique_sort(const std::vector<int>& input)
    {
        std::set<int> result(input.begin(), input.end());

        return std::vector<int>(result.begin(), result.end());
    }

    BootstrappingConfig::BootstrappingConfig(int CtoS, int StoC, int taylor,
                                             bool less_key_mode)
        : CtoS_piece_(CtoS), StoC_piece_(StoC), taylor_number_(taylor),
          less_key_mode_(less_key_mode)
    {
        validate();
    }

    // Validation Function Implementation
    void BootstrappingConfig::validate()
    {
        if (CtoS_piece_ < 2 || CtoS_piece_ > 5)
        {
            throw std::out_of_range("CtoS_piece must be in range [2, 5]");
        }
        if (StoC_piece_ < 2 || StoC_piece_ > 5)
        {
            throw std::out_of_range("StoC_piece must be in range [2, 5]");
        }
        if (taylor_number_ < 6 || taylor_number_ > 15)
        {
            throw std::out_of_range("taylor_number must be in range [6, 15]");
        }
    }

    BootstrappingConfigV2::BootstrappingConfigV2(
        EncodingMatrixConfig stc_config, EvalModConfig eval_mod_config,
        EncodingMatrixConfig cts_config)
        : CtoS_piece_(cts_config.piece_), StoC_piece_(stc_config.piece_),
          stc_config_(stc_config), eval_mod_config_(eval_mod_config),
          cts_config_(cts_config)
    {
    }

    //////////////////////////
    //////////////////////////

    std::vector<Data64>
    calculate_last_q_modinv(const std::vector<Modulus64>& prime_vector,
                            const int Q_prime_size, const int P_size)
    {
        std::vector<Data64> last_q_modinv;
        for (int i = 0; i < P_size; i++)
        {
            for (int j = 0; j < (Q_prime_size - 1) - i; j++)
            {
                // TODO: Change here for BFV as well!!!
                Data64 temp_ = prime_vector[prime_vector.size() - 1 - i].value %
                               prime_vector[j].value;
                last_q_modinv.push_back(
                    OPERATOR64::modinv(temp_, prime_vector[j]));
            }
        }

        return last_q_modinv;
    }

    std::vector<Data64>
    calculate_half(const std::vector<Modulus64>& prime_vector, const int P_size)
    {
        std::vector<Data64> half;
        for (int i = 0; i < P_size; i++)
        {
            half.push_back(prime_vector[prime_vector.size() - 1 - i].value >>
                           1);
        }

        return half;
    }

    std::vector<Data64>
    calculate_half_mod(const std::vector<Modulus64>& prime_vector,
                       const std::vector<Data64>& half, const int Q_prime_size,
                       const int P_size)
    {
        std::vector<Data64> half_mod;
        for (int i = 0; i < P_size; i++)
        {
            for (int j = 0; j < (Q_prime_size - 1) - i; j++)
            {
                half_mod.push_back(half[i] % prime_vector[j].value);
            }
        }

        return half_mod;
    }

    std::vector<Data64>
    calculate_factor(const std::vector<Modulus64>& prime_vector,
                     const int Q_size, const int P_size)
    {
        std::vector<Data64> factor;
        for (int i = 0; i < P_size; i++)
        {
            for (int j = 0; j < Q_size; j++)
            {
                factor.push_back(
                    prime_vector[prime_vector.size() - 1 - i].value %
                    prime_vector[j].value);
            }
        }

        return factor;
    }

    //////////////////////////
    //////////////////////////

    std::vector<Data64> calculate_Mi(const std::vector<Modulus64>& prime_vector,
                                     const int size)
    {
        std::vector<Data64> result_Mi(size * size, 0ULL);
        for (int i = 0; i < size; i++)
        {
            mpz_t result;
            mpz_init(result);
            mpz_set_ui(result, 1);

            for (int j = 0; j < size; j++)
            {
                if (i != j)
                {
                    mpz_mul_ui(result, result, prime_vector[j].value);
                }
            }

            size_t mul_size;
            uint64_t* ptr = reinterpret_cast<uint64_t*>(mpz_export(
                NULL, &mul_size, -1, sizeof(uint64_t), 0, 0, result));

            for (int j = 0; j < mul_size; j++)
            {
                result_Mi[(size * i) + j] = ptr[j];
            }

            mpz_clear(result);
            free(ptr);
        }

        return result_Mi;
    }

    std::vector<Data64>
    calculate_Mi_inv(const std::vector<Modulus64>& prime_vector, const int size)
    {
        std::vector<Data64> result_Mi_inv;
        for (int i = 0; i < size; i++)
        {
            Data64 temp = 1;
            for (int j = 0; j < size; j++)
            {
                if (i != j)
                {
                    Data64 inner_prime =
                        prime_vector[j].value % prime_vector[i].value;
                    temp = OPERATOR64::mult(temp, inner_prime, prime_vector[i]);
                }
            }
            result_Mi_inv.push_back(OPERATOR64::modinv(temp, prime_vector[i]));
        }

        return result_Mi_inv;
    }

    std::vector<Data64> calculate_M(const std::vector<Modulus64>& prime_vector,
                                    const int size)
    {
        std::vector<Data64> result_M(size, 0ULL);

        mpz_t result;
        mpz_init(result);
        mpz_set_ui(result, 1);

        for (int i = 0; i < size; i++)
        {
            mpz_mul_ui(result, result, prime_vector[i].value);
        }

        size_t mul_size;
        uint64_t* ptr = reinterpret_cast<uint64_t*>(
            mpz_export(NULL, &mul_size, -1, sizeof(uint64_t), 0, 0, result));

        for (int j = 0; j < mul_size; j++)
        {
            result_M[j] = ptr[j];
        }

        mpz_clear(result);
        free(ptr);

        return result_M;
    }

    std::vector<Data64>
    calculate_upper_half_threshold(const std::vector<Modulus64>& prime_vector,
                                   const int size)
    {
        std::vector<Data64> result_upper_half_threshold(size, 0ULL);

        mpz_t result;
        mpz_init(result);
        mpz_set_ui(result, 1);

        for (int i = 0; i < size; i++)
        {
            mpz_mul_ui(result, result, prime_vector[i].value);
        }

        mpz_add_ui(result, result, 1);
        mpz_div_2exp(result, result, 1);

        size_t mul_size;
        uint64_t* ptr = reinterpret_cast<uint64_t*>(
            mpz_export(NULL, &mul_size, -1, sizeof(uint64_t), 0, 0, result));

        for (int j = 0; j < mul_size; j++)
        {
            result_upper_half_threshold[j] = ptr[j];
        }

        mpz_clear(result);
        free(ptr);

        return result_upper_half_threshold;
    }

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

} // namespace heongpu