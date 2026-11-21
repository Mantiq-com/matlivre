//! `jbtest` and `dwtest`.

use super::util::*;
use crate::builtins::Ret;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::value::*;

pub const NAMES: &[&str] = &["jbtest", "dwtest", "chi2gof", "vartestn", "ansaribradley", "lillietest"];

const JB_ALPHAS: [f64; 17] = [0.001, 0.0016681, 0.0027826, 0.0046416, 0.0077426, 0.01, 0.012915, 0.021544, 0.025, 0.035938, 0.05, 0.059948, 0.1, 0.15, 0.2, 0.25, 0.50];
const JB_SIZES: [f64; 32] = [
    3.0, 4.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 125.0, 150.0, 175.0, 200.0, 250.0, 300.0, 400.0, 500.0, 800.0, 1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0,
];
const JB_TABLE: [[f64; 17]; 32] = [
    [0.5312, 0.5312, 0.5312, 0.5312, 0.5312, 0.5312, 0.5311, 0.531, 0.5309, 0.5305, 0.5297, 0.529, 0.5251, 0.5176, 0.5074, 0.4946, 0.4063],
    [0.9606, 0.959, 0.9564, 0.952, 0.9448, 0.9396, 0.9329, 0.9133, 0.9056, 0.8817, 0.8519, 0.8316, 0.7553, 0.6721, 0.6303, 0.5947, 0.4739],
    [1.8289, 1.8053, 1.7727, 1.7276, 1.6661, 1.6275, 1.5828, 1.4723, 1.4343, 1.3294, 1.2185, 1.1516, 0.9442, 0.7945, 0.7302, 0.6878, 0.5285],
    [10.9719, 9.843, 8.6751, 7.4815, 6.2927, 5.7077, 5.135, 4.0456, 3.7474, 3.0691, 2.5239, 2.2555, 1.6231, 1.2821, 1.1235, 1.0198, 0.6951],
    [19.5425, 16.7207, 14.0454, 11.5527, 9.2814, 8.2365, 7.2613, 5.5238, 5.0729, 4.0773, 3.2985, 2.9215, 2.0533, 1.5965, 1.3779, 1.2336, 0.7916],
    [25.0722, 21.0001, 17.2981, 13.9762, 11.0602, 9.7531, 8.5489, 6.4401, 5.8998, 4.7176, 3.8011, 3.3596, 2.3505, 1.8239, 1.5631, 1.3885, 0.8577],
    [28.4885, 23.6332, 19.2818, 15.4635, 12.1658, 10.7058, 9.3658, 7.0389, 6.4463, 5.1501, 4.1494, 3.6684, 2.5707, 1.9986, 1.7063, 1.5079, 0.9071],
    [30.6106, 25.2501, 20.5242, 16.4088, 12.8805, 11.3263, 9.9063, 7.4472, 6.8222, 5.4578, 4.4039, 3.8968, 2.7431, 2.139, 1.8216, 1.604, 0.9457],
    [31.9343, 26.2695, 21.3105, 17.0176, 13.3572, 11.7488, 10.2818, 7.7395, 7.0942, 5.6856, 4.5973, 4.0734, 2.8827, 2.2547, 1.9172, 1.6835, 0.9771],
    [32.7514, 26.9025, 21.8016, 17.4145, 13.6749, 12.0355, 10.5399, 7.952, 7.2949, 5.8598, 4.7481, 4.2126, 2.9987, 2.3524, 1.998, 1.7508, 1.0033],
    [33.2273, 27.2686, 22.103, 17.6594, 13.8848, 12.2302, 10.7201, 8.1098, 7.4469, 5.9945, 4.8689, 4.3258, 3.0973, 2.4361, 2.0674, 1.8085, 1.0256],
    [33.4825, 27.4765, 22.2863, 17.8172, 14.0345, 12.3739, 10.8586, 8.234, 7.5661, 6.1045, 4.9697, 4.4214, 3.1834, 2.5097, 2.1283, 1.8592, 1.0449],
    [33.5009, 27.5458, 22.379, 17.9431, 14.1788, 12.5255, 11.0182, 8.3968, 7.7287, 6.262, 5.1203, 4.5688, 3.3246, 2.6316, 2.2297, 1.9434, 1.0768],
    [33.261, 27.3738, 22.2848, 17.9209, 14.2092, 12.5801, 11.0884, 8.4921, 7.8286, 6.3696, 5.2305, 4.6803, 3.4374, 2.7298, 2.3115, 2.0114, 1.1023],
    [32.8742, 27.0975, 22.1135, 17.8305, 14.1911, 12.5841, 11.1155, 8.552, 7.8938, 6.4462, 5.3134, 4.7662, 3.5292, 2.8104, 2.3789, 2.0674, 1.1231],
    [32.3418, 26.7278, 21.877, 17.6977, 14.1283, 12.5542, 11.1112, 8.585, 7.9361, 6.5028, 5.3796, 4.8378, 3.6071, 2.8788, 2.4361, 2.1151, 1.1408],
    [31.8884, 26.399, 21.6407, 17.551, 14.0491, 12.5067, 11.0882, 8.6004, 7.9589, 6.5429, 5.4314, 4.8957, 3.673, 2.9372, 2.4851, 2.1557, 1.1557],
    [30.6089, 25.4748, 21.01, 17.146, 13.8271, 12.3565, 10.9976, 8.602, 7.9816, 6.6071, 5.5277, 5.0071, 3.8025, 3.0523, 2.5819, 2.2363, 1.1852],
    [29.429, 24.6043, 20.3979, 16.7425, 13.5878, 12.1804, 10.8779, 8.5719, 7.9718, 6.6387, 5.5919, 5.0863, 3.8987, 3.1384, 2.6543, 2.2964, 1.2072],
    [28.4225, 23.8582, 19.8658, 16.3846, 13.3664, 12.0139, 10.7589, 8.5307, 7.9496, 6.6571, 5.6408, 5.1474, 3.9732, 3.205, 2.7106, 2.3436, 1.2243],
    [27.514, 23.1914, 19.3922, 16.0635, 13.1643, 11.8602, 10.6494, 8.4904, 7.9259, 6.6688, 5.6783, 5.1957, 4.0327, 3.2584, 2.7559, 2.3812, 1.2382],
    [26.0239, 22.0764, 18.5877, 15.5123, 12.8124, 11.5923, 10.4542, 8.4114, 7.8751, 6.6803, 5.7343, 5.2681, 4.1224, 3.3402, 2.825, 2.4388, 1.2591],
    [24.8353, 21.1758, 17.9384, 15.0633, 12.5206, 11.3653, 10.2839, 8.3384, 7.8276, 6.6842, 5.7728, 5.3194, 4.1873, 3.3988, 2.8749, 2.4808, 1.2744],
    [23.0771, 19.8542, 16.9669, 14.3843, 12.0743, 11.0157, 10.0226, 8.2266, 7.7526, 6.6877, 5.8248, 5.3889, 4.2748, 3.479, 2.9434, 2.5382, 1.2956],
    [21.817, 18.9068, 16.2666, 13.8903, 11.7469, 10.7604, 9.8318, 8.1454, 7.6984, 6.6882, 5.8581, 5.4336, 4.332, 3.5318, 2.9886, 2.576, 1.3097],
    [19.5539, 17.1748, 14.9865, 12.9765, 11.1441, 10.2914, 9.4852, 8.0027, 7.604, 6.6858, 5.9096, 5.5043, 4.4246, 3.618, 3.0631, 2.6388, 1.3334],
    [18.6707, 16.4879, 14.4745, 12.6143, 10.9072, 10.1106, 9.3521, 7.9501, 7.5691, 6.6847, 5.9282, 5.5296, 4.458, 3.65, 3.0907, 2.6624, 1.3421],
    [18.0126, 15.983, 14.099, 12.3493, 10.7362, 9.9803, 9.2587, 7.9127, 7.5441, 6.6832, 5.9408, 5.5471, 4.4814, 3.6718, 3.1099, 2.6787, 1.3484],
    [17.519, 15.6036, 13.8149, 12.1501, 10.6072, 9.8811, 9.1869, 7.8828, 7.5239, 6.68, 5.9482, 5.5587, 4.4979, 3.6876, 3.1237, 2.6904, 1.353],
    [17.1305, 15.3029, 13.5883, 11.9883, 10.5059, 9.8072, 9.1333, 7.8611, 7.5089, 6.6789, 5.9546, 5.5677, 4.5105, 3.6997, 3.1345, 2.6996, 1.3564],
    [16.8158, 15.06, 13.4102, 11.8682, 10.4296, 9.7478, 9.0911, 7.8446, 7.4979, 6.6777, 5.96, 5.5755, 4.5214, 3.7101, 3.1433, 2.7072, 1.3595],
    [16.5585, 14.8633, 13.2657, 11.7666, 10.3636, 9.6981, 9.0549, 7.8296, 7.4872, 6.6763, 5.9635, 5.5806, 4.5289, 3.7173, 3.1499, 2.7129, 1.3619],
];

const LILLIE_Q1: [[f64; 12]; 17] = [
    [0.4328, 0.4308, 0.4292, 0.4217, 0.4131, 0.4065, 0.4008, 0.3754, 0.3453, 0.3213, 0.3028, 0.2581, ],
    [0.4386, 0.4335, 0.4291, 0.4128, 0.3966, 0.3851, 0.3756, 0.3431, 0.3189, 0.3026, 0.2893, 0.2332, ],
    [0.4225, 0.4152, 0.4091, 0.3878, 0.3703, 0.3599, 0.3521, 0.3236, 0.2973, 0.281, 0.2688, 0.2187, ],
    [0.4006, 0.3931, 0.3873, 0.3677, 0.3506, 0.3397, 0.3317, 0.3041, 0.2802, 0.2643, 0.2523, 0.2062, ],
    [0.3828, 0.3749, 0.3692, 0.3494, 0.3326, 0.3223, 0.3146, 0.288, 0.2651, 0.2502, 0.2387, 0.1947, ],
    [0.3656, 0.3582, 0.3524, 0.3332, 0.3171, 0.3071, 0.2997, 0.274, 0.252, 0.2379, 0.2271, 0.1851, ],
    [0.3509, 0.3428, 0.3377, 0.3185, 0.3034, 0.2938, 0.2866, 0.262, 0.241, 0.2274, 0.2171, 0.1767, ],
    [0.3374, 0.3305, 0.3252, 0.3067, 0.2915, 0.2823, 0.2754, 0.2515, 0.2312, 0.2181, 0.2082, 0.1695, ],
    [0.3253, 0.3184, 0.3131, 0.2955, 0.2808, 0.2714, 0.2647, 0.2418, 0.2224, 0.2098, 0.2002, 0.1631, ],
    [0.3147, 0.3074, 0.3021, 0.2849, 0.2706, 0.2618, 0.2553, 0.2333, 0.2145, 0.2025, 0.1932, 0.1574, ],
    [0.3033, 0.297, 0.2923, 0.2756, 0.2619, 0.2536, 0.2473, 0.2257, 0.2075, 0.1958, 0.1868, 0.1522, ],
    [0.296, 0.2893, 0.2844, 0.2673, 0.2539, 0.2457, 0.2398, 0.2189, 0.2012, 0.1898, 0.1811, 0.1475, ],
    [0.2884, 0.2821, 0.2772, 0.2604, 0.2472, 0.239, 0.2331, 0.2126, 0.1953, 0.1843, 0.1759, 0.1433, ],
    [0.28, 0.2736, 0.2688, 0.2532, 0.2403, 0.2324, 0.2266, 0.2068, 0.19, 0.1792, 0.1711, 0.1393, ],
    [0.2732, 0.2667, 0.2621, 0.2465, 0.2341, 0.2265, 0.2208, 0.2013, 0.185, 0.1746, 0.1666, 0.1358, ],
    [0.266, 0.2602, 0.2556, 0.2408, 0.2285, 0.2209, 0.2155, 0.1965, 0.1806, 0.1703, 0.1625, 0.1324, ],
    [0.2603, 0.2545, 0.2503, 0.2352, 0.2232, 0.2159, 0.2105, 0.192, 0.1763, 0.1663, 0.1587, 0.1294, ],
];
const LILLIE_C1: [[f64; 12]; 3] = [
    [1.239, 1.21, 1.189, 1.115, 1.058, 1.023, 0.9964, 0.9092, 0.8358, 0.7893, 0.7539, 0.6183, ],
    [0.1725, 0.1909, 0.1978, 0.158, 0.1582, 0.1639, 0.154, 0.1652, 0.1626, 0.1654, 0.1678, 0.1671, ],
    [0.7237, 0.5838, 0.518, 0.5654, 0.4871, 0.414, 0.418, 0.2731, 0.2215, 0.1729, 0.1323, 0.04738, ],
];
const LILLIE_Q2: [[f64; 12]; 17] = [
    [0.6212, 0.6096, 0.6019, 0.5788, 0.5575, 0.543, 0.5314, 0.4843, 0.4444, 0.42, 0.4009, 0.3165, ],
    [0.5802, 0.5705, 0.5629, 0.5369, 0.5129, 0.4966, 0.4845, 0.4422, 0.4046, 0.3795, 0.3605, 0.2877, ],
    [0.5507, 0.5392, 0.5312, 0.5008, 0.4757, 0.4599, 0.4484, 0.4085, 0.3732, 0.35, 0.332, 0.2645, ],
    [0.519, 0.5068, 0.4988, 0.4702, 0.4466, 0.4317, 0.4202, 0.3813, 0.3485, 0.3266, 0.3099, 0.2458, ],
    [0.4909, 0.4798, 0.4712, 0.4437, 0.4206, 0.4065, 0.3961, 0.3589, 0.3274, 0.307, 0.2913, 0.2309, ],
    [0.468, 0.4563, 0.448, 0.4216, 0.3997, 0.3861, 0.3757, 0.3405, 0.3102, 0.2907, 0.2759, 0.2186, ],
    [0.4475, 0.4368, 0.4289, 0.4025, 0.3805, 0.3673, 0.3577, 0.3242, 0.2955, 0.2768, 0.2626, 0.2082, ],
    [0.4307, 0.4199, 0.4123, 0.3865, 0.3655, 0.3528, 0.3431, 0.3103, 0.2826, 0.2646, 0.251, 0.1991, ],
    [0.4138, 0.4032, 0.3963, 0.3716, 0.3513, 0.3391, 0.3298, 0.2984, 0.2716, 0.2544, 0.2412, 0.1913, ],
    [0.3977, 0.388, 0.3813, 0.3579, 0.3388, 0.3267, 0.3177, 0.287, 0.2614, 0.2449, 0.2323, 0.1842, ],
    [0.3864, 0.3765, 0.3699, 0.3467, 0.3276, 0.3158, 0.3072, 0.2777, 0.2526, 0.2365, 0.2243, 0.1778, ],
    [0.3738, 0.3641, 0.3574, 0.3351, 0.3172, 0.3058, 0.2974, 0.269, 0.2445, 0.229, 0.2172, 0.1721, ],
    [0.3628, 0.354, 0.347, 0.3251, 0.3071, 0.2962, 0.2881, 0.2604, 0.237, 0.2219, 0.2104, 0.1668, ],
    [0.3524, 0.3439, 0.3375, 0.3166, 0.2993, 0.2883, 0.2803, 0.2532, 0.2303, 0.2157, 0.2046, 0.1622, ],
    [0.3441, 0.3357, 0.3295, 0.3079, 0.2911, 0.2806, 0.2729, 0.2466, 0.2243, 0.2101, 0.1991, 0.1577, ],
    [0.3364, 0.3278, 0.3215, 0.3006, 0.284, 0.2738, 0.2661, 0.2404, 0.2186, 0.2046, 0.194, 0.1538, ],
    [0.3288, 0.3197, 0.3138, 0.2938, 0.2774, 0.2672, 0.2597, 0.2344, 0.2132, 0.1997, 0.1893, 0.1501, ],
];
const LILLIE_C2: [[f64; 12]; 3] = [
    [1.544, 1.5, 1.472, 1.377, 1.299, 1.249, 1.214, 1.094, 0.9954, 0.9322, 0.8843, 0.7056, ],
    [0.3558, 0.3058, 0.3233, 0.3535, 0.3127, 0.2727, 0.2549, 0.201, 0.1871, 0.1713, 0.1592, 0.1517, ],
    [-0.112, 0.03913, -0.0726, -0.327, -0.2395, -0.1372, -0.09488, 0.006818, -0.002656, 0.01954, 0.04151, 0.009159, ],
];
const LILLIE_Q3: [[f64; 12]; 17] = [
    [0.4668, 0.4628, 0.4594, 0.4439, 0.4309, 0.4221, 0.415, 0.3847, 0.3502, 0.333, 0.3232, 0.2758, ],
    [0.4497, 0.4418, 0.4363, 0.4153, 0.3961, 0.3844, 0.3765, 0.3512, 0.3273, 0.31, 0.2957, 0.246, ],
    [0.4207, 0.4122, 0.4066, 0.3876, 0.3715, 0.3611, 0.3532, 0.3244, 0.3015, 0.2865, 0.2746, 0.2248, ],
    [0.3984, 0.3913, 0.3851, 0.3648, 0.3483, 0.3384, 0.331, 0.3051, 0.2818, 0.267, 0.2558, 0.2103, ],
    [0.3784, 0.3707, 0.3648, 0.3459, 0.33, 0.32, 0.3127, 0.2875, 0.266, 0.2517, 0.2408, 0.1982, ],
    [0.3625, 0.3548, 0.3487, 0.3296, 0.3142, 0.3045, 0.2974, 0.273, 0.2523, 0.2388, 0.2283, 0.1876, ],
    [0.3473, 0.3398, 0.334, 0.3151, 0.3004, 0.2911, 0.2842, 0.2606, 0.2405, 0.2276, 0.2176, 0.1787, ],
    [0.3312, 0.3244, 0.3192, 0.3021, 0.2875, 0.2786, 0.2721, 0.2496, 0.2303, 0.2178, 0.2083, 0.171, ],
    [0.3193, 0.3127, 0.3079, 0.291, 0.2772, 0.2686, 0.2622, 0.2401, 0.2215, 0.2094, 0.2002, 0.1642, ],
    [0.3102, 0.3034, 0.298, 0.2813, 0.2673, 0.2589, 0.2527, 0.2317, 0.2136, 0.2019, 0.193, 0.1582, ],
    [0.2981, 0.292, 0.2875, 0.2713, 0.258, 0.2501, 0.2442, 0.2237, 0.2063, 0.195, 0.1864, 0.1528, ],
    [0.2896, 0.2835, 0.2789, 0.2634, 0.2506, 0.2426, 0.2367, 0.2167, 0.1998, 0.1888, 0.1804, 0.1479, ],
    [0.2825, 0.2761, 0.271, 0.256, 0.2434, 0.2357, 0.23, 0.2103, 0.1939, 0.1833, 0.1752, 0.1436, ],
    [0.2758, 0.2694, 0.2647, 0.2493, 0.2367, 0.2293, 0.2238, 0.2048, 0.1886, 0.1782, 0.1702, 0.1396, ],
    [0.2683, 0.262, 0.2575, 0.2426, 0.2304, 0.2231, 0.2177, 0.1993, 0.1836, 0.1735, 0.1658, 0.1358, ],
    [0.2609, 0.2553, 0.2511, 0.2365, 0.2248, 0.2177, 0.2125, 0.1944, 0.1791, 0.1693, 0.1617, 0.1324, ],
    [0.2552, 0.2498, 0.2458, 0.2314, 0.2198, 0.2127, 0.2075, 0.1897, 0.1747, 0.1651, 0.1578, 0.1293, ],
];
const LILLIE_C3: [[f64; 12]; 3] = [
    [1.214, 1.186, 1.166, 1.096, 1.039, 1.005, 0.9796, 0.895, 0.8242, 0.7794, 0.7453, 0.6131, ],
    [0.1899, 0.2157, 0.2329, 0.2218, 0.1974, 0.1833, 0.1812, 0.1686, 0.1612, 0.163, 0.1645, 0.1596, ],
    [0.6016, 0.4132, 0.3075, 0.2306, 0.2498, 0.2481, 0.2217, 0.1801, 0.1341, 0.09096, 0.05882, -0.01397, ],
];

/// Shape-preserving cubic interpolation (the slopes of Fritsch and Carlson that MATLAB's `pchip`
/// uses) of `(x, y)` at `q`; `x` increases strictly.
fn pchip_at(x: &[f64], y: &[f64], q: f64) -> f64 {
    let n = x.len();
    let h: Vec<f64> = (0..n - 1).map(|i| x[i + 1] - x[i]).collect();
    let del: Vec<f64> = (0..n - 1).map(|i| (y[i + 1] - y[i]) / h[i]).collect();
    let mut d = vec![0.0; n];
    if n == 2 {
        d[0] = del[0];
        d[1] = del[0];
    } else {
        for k in 1..n - 1 {
            if del[k - 1] * del[k] > 0.0 {
                let hs = h[k - 1] + h[k];
                let w1 = (h[k - 1] + hs) / (3.0 * hs);
                let w2 = (hs + h[k]) / (3.0 * hs);
                let (dmax, dmin) = (del[k - 1].abs().max(del[k].abs()), del[k - 1].abs().min(del[k].abs()));
                d[k] = dmin / (w1 * (del[k - 1] / dmax) + w2 * (del[k] / dmax));
            }
        }
        let end = |h0: f64, h1: f64, d0: f64, d1: f64| {
            let mut s = ((2.0 * h0 + h1) * d0 - h0 * d1) / (h0 + h1);
            if s.signum() != d0.signum() {
                s = 0.0;
            } else if d0.signum() != d1.signum() && s.abs() > 3.0 * d0.abs() {
                s = 3.0 * d0;
            }
            s
        };
        d[0] = end(h[0], h[1], del[0], del[1]);
        d[n - 1] = end(h[n - 2], h[n - 3], del[n - 2], del[n - 3]);
    }
    let mut k = 0;
    while k + 2 < n && q >= x[k + 1] {
        k += 1;
    }
    let t = q - x[k];
    let c = (3.0 * del[k] - 2.0 * d[k] - d[k + 1]) / h[k];
    let b = (d[k] - 2.0 * del[k] + d[k + 1]) / (h[k] * h[k]);
    y[k] + t * (d[k] + t * (c + t * b))
}

/// The critical values of the Jarque-Bera statistic for sample size `n` at the tabulated levels.
fn jb_critical_values(n: f64) -> Vec<f64> {
    // Sizes enter as 1/n; the last row is the limit (chi-square with two degrees of freedom).
    let mut inv: Vec<f64> = vec![0.0];
    inv.extend(JB_SIZES.iter().rev().map(|s| 1.0 / s));
    (0..17)
        .map(|j| {
            let mut col = vec![((-2.0 * JB_ALPHAS[j].ln()) * 1e4).round() / 1e4];
            col.extend(JB_TABLE.iter().rev().map(|row| row[j]));
            pchip_at(&inv, &col, 1.0 / n)
        })
        .collect()
}

fn optional<'a>(a: &'a [Value], i: usize) -> Option<&'a Value> {
    match a.get(i) {
        Some(Value::Num(m)) if m.is_empty() => None,
        other => other,
    }
}

fn jbtest(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    let x = match a.first() {
        Some(Value::Num(m)) => m,
        Some(_) => return rt("Input X must be a vector of real values."),
        None => return rt("Not enough input arguments."),
    };
    let is_vec = x.dims.len() == 2 && (x.dims[0] == 1 || x.dims[1] == 1);
    if !is_vec {
        return rt("Input X must be a vector of real values.");
    }
    let xs: Vec<f64> = x.data.iter().copied().filter(|v| !v.is_nan()).collect();
    let n = xs.len();
    if n < 2 {
        return rt("Sample vector X must have at least 2 valid observations.");
    }
    if n == 2 {
        return Ok(Some(vec![Value::num(0.0), Value::num(1.0), Value::num(1.0 / 3.0), Value::num(f64::NAN)]));
    }
    let alpha = match optional(a, 1) {
        Some(Value::Num(m)) if m.data.len() == 1 && 0.0 < m.data[0] && m.data[0] < 1.0 => m.data[0],
        Some(_) => return rt("Significance level ALPHA must be a scalar between 0 and 1."),
        None => 0.05,
    };
    let mctol = match optional(a, 2) {
        Some(Value::Num(m)) if m.data.len() == 1 && m.data[0] > 0.0 => Some(m.data[0]),
        Some(_) => return rt("Monte-Carlo standard error tolerance MCTOL must be a positive scalar value."),
        None => None,
    };
    let nf = n as f64;
    let mean = xs.iter().sum::<f64>() / nf;
    let sd1 = (xs.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / nf).sqrt();
    let z: Vec<f64> = xs.iter().map(|v| (v - mean) / sd1).collect();
    let skew = z.iter().map(|v| v.powi(3)).sum::<f64>() / nf;
    let kurt = z.iter().map(|v| v.powi(4)).sum::<f64>() / nf - 3.0;
    let jb = nf * (skew * skew / 6.0 + kurt * kurt / 24.0);
    let (cv, p);
    if let Some(tol) = mctol {
        let (c, pp) = jb_monte_carlo(it, jb, n, alpha, tol)?;
        cv = c;
        p = pp;
    } else {
        if alpha < JB_ALPHAS[0] || JB_ALPHAS[16] < alpha {
            return rt("Significance level ALPHA is outside the range of the tabulated values. Use a value in the interval [0.001, 0.5], or use the MCTOL input argument.");
        }
        let cvs = jb_critical_values(nf);
        cv = pchip_at(&JB_ALPHAS, &cvs, alpha);
        let mut pv = f64::NAN;
        if nargout > 1 {
            if jb < cvs[16] {
                it.emit_warning("P is greater than the largest tabulated value, returning 0.5.".into(), "stats:jbtest:PTooBig".into());
                pv = JB_ALPHAS[16];
            } else if cvs[0] < jb {
                it.emit_warning("P is less than the smallest tabulated value, returning 0.001.".into(), "stats:jbtest:PTooSmall".into());
                pv = JB_ALPHAS[0];
            } else if jb == cvs[0] {
                pv = JB_ALPHAS[0];
            } else if jb.is_nan() {
                pv = 0.0;
            } else {
                let i = cvs.iter().position(|&c| jb > c).unwrap();
                let (mut lo, mut hi) = (JB_ALPHAS[i - 1], JB_ALPHAS[i]);
                let f = |t: f64| pchip_at(&JB_ALPHAS, &cvs, t) - jb;
                let flo = f(lo);
                for _ in 0..200 {
                    let mid = 0.5 * (lo + hi);
                    if (f(mid) > 0.0) == (flo > 0.0) {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                pv = 0.5 * (lo + hi);
            }
        }
        p = pv;
    }
    let h = if jb > cv { 1.0 } else { 0.0 } + if jb.is_nan() { 1.0 } else { 0.0 };
    let single = is_single(x);
    let sv = |v: f64| Value::Num(finish(Matrix::scalar(v), single));
    Ok(Some(vec![Value::num(h), sv(p), sv(jb), Value::num(cv)]))
}

fn jb_monte_carlo(it: &mut Interpreter, jb: f64, n: usize, alpha: f64, tol: f64) -> Res<(f64, f64)> {
    let vartol = tol * tol;
    let (mut crit, mut p, mut reps_tot) = (0.0, 0.0, 0.0f64);
    let mut reps_min = 1000.0f64;
    loop {
        let old = reps_tot;
        let reps = (reps_min - old).ceil().max(1.0) as usize;
        let mut stats = Vec::with_capacity(reps);
        for _ in 0..reps {
            let r = it.call_function("randn", vec![Value::num(1.0), Value::num(n as f64)], 1)?;
            let Value::Num(m) = &r[0] else { return rt("randn") };
            let nf = n as f64;
            let mean = m.data.iter().sum::<f64>() / nf;
            let s2 = m.data.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / nf;
            let skew = m.data.iter().map(|v| (v - mean).powi(3)).sum::<f64>() / (s2.powf(1.5) * nf);
            let kurt = m.data.iter().map(|v| (v - mean).powi(4)).sum::<f64>() / (s2 * s2 * nf) - 3.0;
            stats.push(nf * (skew * skew / 6.0 + kurt * kurt / 24.0));
        }
        let q = it.call_function("prctile", vec![Value::Num(Matrix::from_dims(vec![reps, 1], stats.clone())), Value::num(100.0 * (1.0 - alpha))], 1)?;
        let crit_mc = q[0].scalar()?;
        let p_mc = stats.iter().filter(|&&s| s > jb).count() as f64 / reps as f64;
        reps_tot = old + reps as f64;
        crit = (old * crit + reps as f64 * crit_mc) / reps_tot;
        p = (old * p + reps as f64 * p_mc) / reps_tot;
        let sepsq = (p * (1.0 - p) / reps_tot).max(1.0 / (reps_tot * reps_tot));
        if sepsq < vartol {
            break;
        }
        reps_min = 1.2 * (reps_tot * sepsq) / vartol;
    }
    if jb.is_nan() {
        p = 0.0;
    }
    Ok((crit, p))
}

// ---------------------------------------------------------------------------------------------
// dwtest

type Mat = Vec<Vec<f64>>;

fn matmul(a: &Mat, b: &Mat) -> Mat {
    let (n, k, m) = (a.len(), b.len(), b[0].len());
    let mut c = vec![vec![0.0; m]; n];
    for i in 0..n {
        for l in 0..k {
            let v = a[i][l];
            for j in 0..m {
                c[i][j] += v * b[l][j];
            }
        }
    }
    c
}

fn transpose(a: &Mat) -> Mat {
    (0..a[0].len()).map(|j| a.iter().map(|r| r[j]).collect()).collect()
}

fn trace(a: &Mat) -> f64 {
    (0..a.len()).map(|i| a[i][i]).sum()
}

fn inverse(a: &Mat) -> Mat {
    let n = a.len();
    let mut m: Mat = (0..n).map(|i| a[i].iter().copied().chain((0..n).map(|j| if i == j { 1.0 } else { 0.0 })).collect()).collect();
    for c in 0..n {
        let piv = (c..n).max_by(|&i, &j| m[i][c].abs().partial_cmp(&m[j][c].abs()).unwrap()).unwrap();
        m.swap(c, piv);
        let d = m[c][c];
        for v in m[c].iter_mut() {
            *v /= d;
        }
        for r in 0..n {
            if r != c {
                let f = m[r][c];
                if f != 0.0 {
                    for j in 0..2 * n {
                        let t = m[c][j];
                        m[r][j] -= f * t;
                    }
                }
            }
        }
    }
    m.into_iter().map(|r| r[n..].to_vec()).collect()
}

/// P(DW < dw) from the eigenvalues `lam` (ascending) of the quadratic form, by numerical
/// integration of the characteristic function (Pan's method as described by Farebrother).
fn pan_cdf(lam: &[f64], dw: f64, n: usize) -> f64 {
    let m = lam.len() as i64;
    let Some(first) = lam.iter().position(|&l| l >= dw) else { return 0.0 };
    let mut mu = first as i64;
    let mut h = m - mu;
    let (d, mut j1, j2, mut j3, mut j4, k);
    if h > mu {
        d = 2;
        h = mu;
        k = -1.0;
        j1 = 0;
        j2 = 2;
        j3 = 3;
        j4 = 1;
    } else {
        d = -2;
        mu += 1;
        j1 = m - 2;
        j2 = m - 1;
        j3 = m + 1;
        j4 = m;
        k = 1.0;
    }
    let step = std::f64::consts::FRAC_PI_2 / n as f64;
    let mut sum = (k + 1.0) / 2.0;
    let mut sgn = k / n as f64;
    let at = |i: i64| lam[(i - 1) as usize];
    let mut start = j2;
    let mut l1 = h - 2 * (h / 2);
    while l1 >= 0 {
        let mut l2 = start;
        while (d > 0 && l2 <= mu) || (d < 0 && l2 >= mu) {
            let top = at(j4);
            let other = if l2 == 0 { dw } else { at(l2) };
            let u = 0.5 * (top + other);
            let v = 0.5 * (top - other);
            let mut acc = 0.0;
            let mut i = 1;
            while i <= 2 * n - 1 {
                let y = u - v * (i as f64 * step).cos();
                let mut prod = 1.0;
                for t in (1..=j1).chain(j3..=m) {
                    prod *= (y - dw) / (y - at(t));
                }
                acc += prod.abs().sqrt();
                i += 2;
            }
            sgn = -sgn;
            sum += sgn * acc;
            j1 += d;
            j3 += d;
            j4 += d;
            l2 += d;
        }
        if d == 2 {
            j3 -= 1;
        } else {
            j1 += 1;
        }
        start = 0;
        mu = 0;
        l1 -= 1;
    }
    sum
}

fn dw_approx(it: &mut Interpreter, x: &Mat, dw: f64) -> Res<f64> {
    let (n, p) = (x.len(), x[0].len());
    let xt = transpose(x);
    let q1 = inverse(&matmul(&xt, x));
    let mut b: Mat = Vec::with_capacity(n);
    b.push((0..p).map(|j| x[0][j] - x[1][j]).collect());
    for r in 1..n - 1 {
        b.push((0..p).map(|j| 2.0 * x[r][j] - x[r - 1][j] - x[r + 1][j]).collect());
    }
    b.push((0..p).map(|j| -x[n - 2][j] + x[n - 1][j]).collect());
    let c = matmul(&matmul(&xt, &b), &q1);
    let nu1 = 2.0 * (n as f64 - 1.0) - trace(&c);
    let btb = matmul(&matmul(&transpose(&b), &b), &q1);
    let nu2 = 2.0 * (3.0 * n as f64 - 4.0) - 2.0 * trace(&btb) + trace(&matmul(&c, &c));
    let df = n as f64 - p as f64;
    let mu = nu1 / df;
    let sigma = (2.0 / (df * (df + 2.0)) * (nu2 - nu1 * mu)).sqrt();
    let r = it.call_function("normcdf", vec![Value::num(dw), Value::num(mu), Value::num(sigma)], 1)?;
    r[0].scalar()
}

fn pick(given: &str, choices: &[&str], what: &str) -> Res<String> {
    let g = given.to_lowercase();
    if !g.is_empty() {
        let hits: Vec<&&str> = choices.iter().filter(|c| c.starts_with(&g)).collect();
        if hits.len() == 1 {
            return Ok(hits[0].to_string());
        }
    }
    let list: Vec<String> = choices.iter().map(|c| format!("'{c}'")).collect();
    rt(format!("'{given}' is not a valid value for the '{what}' argument.  Valid values are: {}.", list.join(", ")))
}

fn text_of(v: &Value) -> String {
    if v.is_text() {
        v.text().unwrap_or_default()
    } else {
        String::new()
    }
}

fn dwtest(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    let _ = nargout;
    let (Some(Value::Num(r)), Some(Value::Num(xm))) = (a.first(), a.get(1)) else { return rt("Not enough input arguments.") };
    let (n, p) = (xm.dims[0], xm.dims[1]);
    let mut option = String::new();
    let mut alt = String::new();
    if a.len() > 2 {
        let v = &a[2];
        if !v.is_text() {
            return rt(format!("Input A of class {} and input B of class cell must be cell arrays of character vectors, unless one is a character vector.", v.class_name()));
        }
        let low = text_of(v).to_lowercase();
        if low == "exact" || low == "approximate" {
            option = low;
            if a.len() > 3 {
                alt = text_of(&a[3]);
            }
        } else if a.len() == 3 {
            return rt("The METHOD argument must be 'approximate' or 'exact'.");
        } else {
            let mut i = 2;
            while i < a.len() {
                let key = text_of(&a[i]).to_lowercase();
                let val = a.get(i + 1).map(text_of).unwrap_or_default();
                if !key.is_empty() && "method".starts_with(&key) {
                    option = val;
                } else if !key.is_empty() && "tail".starts_with(&key) {
                    alt = val;
                } else {
                    return rt(format!("Invalid parameter name: {}.", text_of(&a[i])));
                }
                i += 2;
            }
        }
    }
    let option = if option.is_empty() {
        if n < 400 { "exact".to_string() } else { "approximate".to_string() }
    } else {
        pick(&option, &["exact", "approximate"], "method")?
    };
    let alt = if alt.is_empty() { "both".to_string() } else { pick(&alt, &["left", "both", "right"], "tail")? };
    let d: Vec<f64> = r.data.windows(2).map(|w| w[1] - w[0]).collect();
    let dw = d.iter().map(|v| v * v).sum::<f64>() / r.data.iter().map(|v| v * v).sum::<f64>();
    let xrows: Mat = (0..n).map(|i| (0..p).map(|j| xm.data[j * n + i]).collect()).collect();
    let mut pdw = f64::NAN;
    let mut use_exact = option == "exact";
    if use_exact {
        let q = it.call_function("qr", vec![Value::Num(xm.clone())], 2)?;
        let Value::Num(qm) = &q[0] else { return rt("qr") };
        // G = F * Q(:, p+1:n), F the (n-1)-by-n difference matrix.
        let cols = n - p;
        let mut g = vec![0.0; (n - 1) * cols];
        for c in 0..cols {
            for i in 0..n - 1 {
                g[c * (n - 1) + i] = qm.data[(p + c) * n + i + 1] - qm.data[(p + c) * n + i];
            }
        }
        let sv = it.call_function("svd", vec![Value::Num(Matrix::from_dims(vec![n - 1, cols], g))], 1)?;
        let Value::Num(s) = &sv[0] else { return rt("svd") };
        let mut e: Vec<f64> = s.data.iter().map(|v| v * v).collect();
        e.sort_by(|a, b| a.partial_cmp(b).unwrap());
        pdw = pan_cdf(&e, dw, n);
        if pdw > 1.0 || pdw < 0.0 {
            it.emit_warning("Exact p-value is unavailable, using approximate p-value.".into(), "stats:pvaluedw:ExactUnavailable".into());
            use_exact = false;
        }
    }
    if !use_exact {
        pdw = dw_approx(it, &xrows, dw)?;
    }
    let pv = match alt.as_str() {
        "both" => 2.0 * pdw.min(1.0 - pdw),
        "left" => 1.0 - pdw,
        _ => pdw,
    };
    Ok(Some(vec![Value::num(pv), Value::num(dw)]))
}

pub fn call(it: &mut Interpreter, name: &str, a: &[Value], nargout: usize) -> Ret {
    match name {
        "jbtest" => jbtest(it, a, nargout),
        "chi2gof" => chi2gof(it, a, nargout),
        "vartestn" => vartestn(it, a, nargout),
        "ansaribradley" => ansaribradley(it, a, nargout),
        "lillietest" => lillietest(it, a, nargout),
        _ => dwtest(it, a, nargout),
    }
}

// ---------------------------------------------------------------------------------------------
// chi2gof

fn is_empty_val(v: &Value) -> bool {
    matches!(v, Value::Num(m) if m.is_empty())
}

fn row(v: Vec<f64>) -> Value {
    let n = v.len();
    Value::Num(Matrix::from_dims(vec![1, n], v))
}

fn is_vector_dims(d: &[usize]) -> bool {
    d.len() == 2 && (d[0] == 1 || d[1] == 1)
}

/// The bin counts and edges of `x` for a bin specification (nbins, centers or edges).
fn get_bins(it: &mut Interpreter, x: &[f64], freq: &[f64], kind: &str, spec: &Matrix) -> Res<(Vec<f64>, Vec<f64>)> {
    let lo = x.iter().cloned().fold(f64::INFINITY, f64::min);
    let hi = x.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let mut edges: Vec<f64>;
    let use_centers = kind != "edges";
    match kind {
        "nbins" => {
            let nb = spec.data[0];
            let (mut lo, mut hi) = if x.is_empty() { (0.0, 1.0) } else { (lo, hi) };
            if lo == hi {
                lo = lo - (nb / 2.0).floor() - 0.5;
                hi = hi + (nb / 2.0).ceil() - 0.5;
            }
            let bw = (hi - lo) / nb;
            let count = if nb < 0.0 { 0 } else { nb.floor() as usize + 1 };
            edges = (0..count).map(|k| lo + bw * k as f64).collect();
            if edges.is_empty() {
                return rt("Array indices must be positive integers or logical values.");
            }
            if count == 1 {
                return rt("First input IND and third input SZ must satisfy ALL(MAX(IND)<=SZ).");
            }
            let last = edges.len() - 1;
            edges[last] = hi;
        }
        "ctrs" => {
            let c = &spec.data;
            let mut bw: Vec<f64> = c.windows(2).map(|w| w[1] - w[0]).collect();
            if bw.iter().any(|&b| b <= 0.0) {
                return rt("The vector C of bin centers must be increasing.");
            }
            let l = *bw.last().unwrap_or(&f64::NAN);
            bw.push(l);
            edges = vec![c[0] - bw[0] / 2.0];
            edges.extend(c.iter().zip(&bw).map(|(ci, b)| ci + b / 2.0));
            if !x.is_empty() {
                edges[0] = edges[0].min(lo);
                let last = edges.len() - 1;
                edges[last] = edges[last].max(hi);
            }
        }
        _ => {
            edges = spec.data.clone();
            if edges.windows(2).any(|w| w[1] - w[0] <= 0.0) {
                return rt("The vector E of bin edges must be increasing.");
            }
            if !x.is_empty() && (lo < edges[0] || edges[edges.len() - 1] < hi) {
                it.emit_warning("The vector E of bin edges does not span the range of X.".into(), "stats:ecdfhist:DataNotSpanned".into());
            }
        }
    }
    let nb = edges.len() - 1;
    let mut counts = vec![0.0; nb];
    if use_centers {
        for e in edges.iter_mut() {
            *e += eps_of(*e);
        }
        for (v, w) in x.iter().zip(freq) {
            let k = edges[1..nb].iter().filter(|&&e| e <= *v).count();
            counts[k] += w;
        }
    } else {
        for (v, w) in x.iter().zip(freq) {
            if *v < edges[0] || *v > edges[nb] {
                continue;
            }
            let mut k = edges.iter().filter(|&&e| e <= *v).count().saturating_sub(1);
            if k >= nb {
                k = nb - 1;
            }
            counts[k] += w;
        }
    }
    Ok((counts, edges))
}

fn chi2gof(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    let Some(Value::Num(xm)) = a.first() else { return rt("Not enough input arguments.") };
    if !is_vector_dims(&xm.dims) {
        return rt("X must be a vector of real values.");
    }
    let names = ["nbins", "ctrs", "edges", "cdf", "expected", "nparams", "emin", "frequency", "alpha"];
    let mut vals: Vec<Option<Value>> = vec![None; names.len()];
    if (a.len() - 1) % 2 != 0 {
        return rt("Wrong number of arguments.");
    }
    for pair in a[1..].chunks(2) {
        let key = text_of(&pair[0]).to_lowercase();
        let hits: Vec<usize> = names.iter().enumerate().filter(|(_, n)| !key.is_empty() && n.starts_with(&key)).map(|(i, _)| i).collect();
        let exact = names.iter().position(|n| *n == key);
        let idx = match (exact, hits.len()) {
            (Some(i), _) => i,
            (None, 1) => hits[0],
            (None, n) if n > 1 => return rt(format!("Ambiguous parameter name: {}.", text_of(&pair[0]))),
            _ => return rt(format!("Invalid parameter name: {}.", text_of(&pair[0]))),
        };
        vals[idx] = Some(pair[1].clone());
    }
    let get = |i: usize| -> Option<&Value> { vals[i].as_ref().filter(|v| !is_empty_val(v)) };
    let (nbins_v, ctrs_v, edges_v, cdf_v, exp_v, np_v, emin_v, freq_v, alpha_v) = (get(0), get(1), get(2), get(3), get(4), get(5), vals[6].as_ref(), get(7), vals[8].as_ref());
    if [nbins_v, ctrs_v, edges_v].iter().filter(|v| v.is_some()).count() > 1 {
        return rt("You cannot specify more than one of 'nbins', 'ctrs', and 'edges'.");
    }
    if cdf_v.is_some() && exp_v.is_some() {
        return rt("You cannot specify both 'cdf' and 'expected'.");
    }
    let freq_in: Option<Matrix> = match freq_v {
        Some(Value::Num(m)) => Some(m.clone()),
        _ => None,
    };
    if let Some(f) = &freq_in {
        if !is_vector_dims(&f.dims) || f.data.len() != xm.data.len() {
            return rt("The 'frequency' vector must have the same size as X.");
        }
        if f.data.iter().any(|&v| v < 0.0) {
            return rt("The 'frequency' values must be zero or positive.");
        }
    }
    let emin = match emin_v {
        None => 5.0,
        Some(Value::Num(m)) if m.data.len() == 1 && m.data[0] >= 0.0 && m.data[0] == m.data[0].round() => m.data[0],
        _ => return rt("The 'emin' value must be a non-negative integer."),
    };
    let mut nparams: Option<f64> = match np_v {
        None => None,
        Some(Value::Num(m)) if m.data.len() == 1 && m.data[0] >= 0.0 && m.data[0] == m.data[0].round() => Some(m.data[0]),
        _ => return rt("The 'nparams' value must be a non-negative integer."),
    };
    let alpha = match alpha_v {
        None => 0.05,
        Some(Value::Num(m)) if m.data.len() == 1 && m.data[0] > 0.0 && m.data[0] < 1.0 => m.data[0],
        _ => return rt("The 'alpha' value must be a scalar between 0 and 1."),
    };
    let freq0 = freq_in.map(|f| f.data).unwrap_or_else(|| vec![1.0; xm.data.len()]);
    let keep: Vec<bool> = xm.data.iter().zip(&freq0).map(|(a, b)| !(a.is_nan() || b.is_nan())).collect();
    let x: Vec<f64> = xm.data.iter().zip(&keep).filter(|(_, k)| **k).map(|(v, _)| *v).collect();
    let freq: Vec<f64> = freq0.iter().zip(&keep).filter(|(_, k)| **k).map(|(v, _)| *v).collect();
    let (obs, mut edges) = if let Some(Value::Num(c)) = ctrs_v {
        get_bins(it, &x, &freq, "ctrs", c)?
    } else if let Some(Value::Num(e)) = edges_v {
        get_bins(it, &x, &freq, "edges", e)?
    } else {
        let nb = match (nbins_v, exp_v) {
            (Some(Value::Num(m)), _) => m.data[0],
            (_, Some(Value::Num(e))) => e.dims.iter().copied().max().unwrap_or(0) as f64,
            _ => 10.0,
        };
        get_bins(it, &x, &freq, "nbins", &Matrix::scalar(nb))?
    };
    let mut obs = obs;
    let mut nb = obs.len();
    let mut exp: Vec<f64>;
    if let Some(ev) = exp_v {
        let bad = || rt::<()>(format!("The 'expected' value must be a vector with {nb} elements."));
        let Value::Num(e) = ev else { return bad().map(|_| None) };
        if !is_vector_dims(&e.dims) || e.data.len() != nb {
            return bad().map(|_| None);
        }
        if e.data.iter().any(|&v| v < 0.0) {
            return rt("The 'expected' values must be non-negative.");
        }
        exp = e.data.clone();
    } else {
        let sumf: f64 = freq.iter().sum();
        let (spec, mut cargs): (Value, Vec<Value>) = match cdf_v {
            None => {
                let mu = x.iter().zip(&freq).map(|(v, w)| v * w).sum::<f64>() / sumf;
                let var = x.iter().zip(&freq).map(|(v, w)| w * (v - mu) * (v - mu)).sum::<f64>() / sumf;
                if nparams.is_none() {
                    nparams = Some(2.0);
                }
                (Value::chars("normcdf"), vec![Value::num(mu), Value::num((sumf / (sumf - 1.0) * var).sqrt())])
            }
            Some(Value::Cell(c)) => {
                let rest: Vec<Value> = c.data[1..].to_vec();
                if nparams.is_none() {
                    nparams = Some(rest.len() as f64);
                }
                (c.data[0].clone(), rest)
            }
            Some(v) => (v.clone(), vec![]),
        };
        if !(spec.is_text() || matches!(spec, Value::Func(_))) {
            return rt("The 'cdf' value must be a ProbDist object, a function handle, or a cell array containing a function handle.");
        }
        let fname = crate::builtins_funfun::func_text(it, &spec);
        let interior: Vec<f64> = if edges.len() >= 2 { edges[1..edges.len() - 1].to_vec() } else { vec![] };
        let mut args = vec![row(interior)];
        args.append(&mut cargs);
        let res = match &spec {
            Value::Func(f) => it.call_value(f, args, 1),
            other => it.call_function(&other.text().unwrap_or_default(), args, 1),
        };
        let fc = match res {
            Ok(r) => r,
            Err(_) => return rt(format!("Error evaluating cdf function '{fname}'.")),
        };
        let fm = match fc.first() {
            Some(Value::Num(m)) if is_vector_dims(&m.dims) && m.data.len() == nb - 1 => m.clone(),
            _ => return rt(format!("The cdf function '{fname}' returned an incorrect number of values.")),
        };
        let total: f64 = obs.iter().sum();
        let mut prev = 0.0;
        exp = Vec::with_capacity(nb);
        for v in fm.data.iter().chain(std::iter::once(&1.0)) {
            exp.push(total * (v - prev));
            prev = *v;
        }
    }
    if exp.iter().any(|&e| e < emin) {
        let (mut i, mut j) = (0usize, exp.len() - 1);
        while i + 1 < j && (exp[i] < emin || exp[i + 1] < emin || exp[j] < emin || exp[j - 1] < emin) {
            if exp[i] < exp[j] {
                exp[i + 1] += exp[i];
                obs[i + 1] += obs[i];
                i += 1;
            } else {
                exp[j - 1] += exp[j];
                obs[j - 1] += obs[j];
                j -= 1;
            }
        }
        exp = exp[i..=j].to_vec();
        obs = obs[i..=j].to_vec();
        edges.drain(j + 1..edges.len() - 1);
        edges.drain(1..=i);
        nb = exp.len();
        if exp.iter().any(|&e| e < emin) {
            it.emit_warning("After pooling, some bins still have low expected counts.\nThe chi-square approximation may not be accurate".into(), "stats:chi2gof:LowCounts".into());
        }
    }
    let cstat: f64 = obs.iter().zip(&exp).map(|(o, e)| (o - e) * (o - e) / e).sum();
    let df = nb as f64 - 1.0 - nparams.unwrap_or(0.0);
    let (df, p) = if df > 0.0 {
        let r = it.call_function("chi2cdf", vec![Value::num(cstat), Value::num(df), Value::chars("upper")], 1)?;
        (df, r[0].scalar()?)
    } else {
        (0.0, f64::NAN)
    };
    let h = if p <= alpha { 1.0 } else { 0.0 };
    let single = is_single(xm);
    let sv = |v: f64| Value::Num(finish(Matrix::scalar(v), single));
    let srow = |v: Vec<f64>| {
        let n = v.len();
        Value::Num(finish(Matrix::from_dims(vec![1, n], v), single))
    };
    let mut out = vec![sv(h), sv(p)];
    if nargout > 2 {
        let mut st = StructArray::scalar();
        st.set("chi2stat", sv(cstat));
        st.set("df", Value::num(df));
        st.set("edges", row(edges));
        st.set("O", srow(obs));
        st.set("E", srow(exp));
        out.push(Value::Struct(st));
    }
    Ok(Some(out))
}

// ---------------------------------------------------------------------------------------------
// vartestn

/// MATLAB's matching of a value against a list of names: exact (any case) first, then a unique
/// prefix. `what` is used as in "the {what} argument"; with `list` the valid values are named.
fn param_val(v: &Value, choices: &[&str], what: &str, list: bool) -> Res<usize> {
    if !v.is_text() {
        return rt(format!("Invalid '{}' data type. Argument must be a character vector or string scalar.", what.trim_matches('\'')));
    }
    let given = v.text().unwrap_or_default();
    let t = given.to_lowercase();
    if let Some(i) = choices.iter().position(|c| c.to_lowercase() == t) {
        return Ok(i);
    }
    let hits: Vec<usize> = choices.iter().enumerate().filter(|(_, c)| !t.is_empty() && c.to_lowercase().starts_with(&t)).map(|(i, _)| i).collect();
    match hits.len() {
        1 => Ok(hits[0]),
        0 if list => {
            let names: Vec<String> = choices.iter().map(|c| format!("'{c}'")).collect();
            rt(format!("'{given}' is not a valid value for the {what} argument.  Valid values are: {}.", names.join(", ")))
        }
        0 => rt(format!("'{given}' is not a valid value for the {what} argument.")),
        _ => rt(format!("The value '{given}' is a partial match for more than one valid value for the {what} argument.")),
    }
}

/// One-way ANOVA of the values by group label, through `anova1`: (F as anova1 reports it, df1, df2, p).
fn oneway(it: &mut Interpreter, x: &[f64], g: &[usize]) -> Res<(Value, f64, f64, f64)> {
    let n = x.len();
    let xv = Value::Num(Matrix::from_dims(vec![n, 1], x.to_vec()));
    let gv = Value::Num(Matrix::from_dims(vec![n, 1], g.iter().map(|&k| (k + 1) as f64).collect()));
    let r = it.call_function("anova1", vec![xv, gv, Value::chars("off")], 2)?;
    let pv = r[0].scalar()?;
    let Value::Cell(c) = &r[1] else { return rt("anova1") };
    let rows = c.dims[0];
    let cell = |i: usize, j: usize| c.data[(j - 1) * rows + (i - 1)].clone();
    Ok((cell(2, 5), cell(2, 3).scalar()?, cell(3, 3).scalar()?, pv))
}

fn vartestn(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    if a.is_empty() {
        return rt("Not enough input arguments.");
    }
    let mut group: Value = a.get(1).cloned().unwrap_or_else(|| Value::Num(Matrix::from_dims(vec![0, 0], vec![])));
    let mut rest: Vec<Value> = a.get(2..).map(|r| r.to_vec()).unwrap_or_default();
    let mut displayopt: Option<Value> = None;
    let mut testtype: Option<Value> = None;
    if a.len() >= 3 {
        if group.is_text() && group.text().map(|t| !t.is_empty()).unwrap_or(false) {
            rest.insert(0, group.clone());
            group = Value::Num(Matrix::from_dims(vec![0, 0], vec![]));
        }
        let first = &rest[0];
        let ft = if first.is_text() { first.text().unwrap_or_default() } else { String::new() };
        if ft == "on" || ft == "off" || is_empty_val(first) {
            displayopt = Some(first.clone());
            if a.len() >= 4 && rest.len() >= 2 {
                testtype = Some(rest[1].clone());
            }
        } else {
            if rest.len() % 2 != 0 {
                return rt("Wrong number of arguments.");
            }
            for pair in rest.chunks(2) {
                let key = text_of(&pair[0]).to_lowercase();
                if !key.is_empty() && "displayopt".starts_with(&key) {
                    displayopt = Some(pair[1].clone());
                } else if !key.is_empty() && "testtype".starts_with(&key) {
                    testtype = Some(pair[1].clone());
                } else {
                    return rt(format!("Invalid parameter name: {}.", text_of(&pair[0])));
                }
            }
        }
    }
    if let Some(d) = &displayopt {
        if !is_empty_val(d) && !(d.is_text() && matches!(d.text().as_deref(), Ok("on") | Ok("off"))) {
            return rt("DISPLAYOPT must be 'on' or 'off'.");
        }
    }
    let dorobust: usize = match &testtype {
        Some(t) if !is_empty_val(t) => {
            let ttypes = ["robust", "classical", "LeveneAbsolute", "OBrien", "BrownForsythe", "LeveneQuadratic", "Bartlett"];
            param_val(t, &ttypes, "'testtype'", false)? + 1
        }
        _ => 0,
    };
    let Value::Num(xm) = &a[0] else { return rt("X must not be empty") };
    if xm.is_empty() {
        return rt("X must not be empty");
    }
    let (n, m) = (xm.dims[0], xm.data.len() / xm.dims[0].max(1));
    let is_col = xm.dims.len() == 2 && xm.dims[1] == 1;
    let (xs, gv): (Vec<f64>, Value) = if is_empty_val(&group) {
        let g: Vec<f64> = (0..m).flat_map(|j| std::iter::repeat((j + 1) as f64).take(n)).collect();
        (xm.data.clone(), Value::Num(Matrix::from_dims(vec![n * m, 1], g)))
    } else if is_col {
        let glen = match &group {
            Value::Num(g) => g.data.len(),
            Value::Cell(c) => c.data.len(),
            other if other.is_text() => other.text().map(|t| t.chars().count()).unwrap_or(0),
            _ => 0,
        };
        let gdims = match &group {
            Value::Cell(c) => c.dims.clone(),
            Value::Num(g) => g.dims.clone(),
            _ => vec![glen, 1],
        };
        // A character row is a single label; a char matrix has one label per row.
        let _ = gdims;
        if glen != n && !(group.is_text() && group.text().map(|_| false).unwrap_or(false)) {
            return rt("X and GROUP must have the same length.");
        }
        (xm.data.clone(), group.clone())
    } else {
        return rt("Cannot specify a GROUP argument unless X is a column vector.");
    };
    // Drop NaN values and missing group labels.
    let keep: Vec<bool> = (0..xs.len())
        .map(|i| {
            let gmiss = match &gv {
                Value::Num(g) => g.data.get(i).map(|v| v.is_nan()).unwrap_or(false),
                Value::Cell(c) => c.data.get(i).map(|v| v.text().map(|t| t.is_empty()).unwrap_or(false)).unwrap_or(false),
                _ => false,
            };
            !(xs[i].is_nan() || gmiss)
        })
        .collect();
    let x: Vec<f64> = xs.iter().zip(&keep).filter(|(_, k)| **k).map(|(v, _)| *v).collect();
    let gsel = match &gv {
        Value::Num(g) => Value::Num(Matrix::from_dims(vec![x.len(), 1], g.data.iter().zip(&keep).filter(|(_, k)| **k).map(|(v, _)| *v).collect())),
        Value::Cell(c) => {
            let d: Vec<Value> = c.data.iter().zip(&keep).filter(|(_, k)| **k).map(|(v, _)| v.clone()).collect();
            let l = d.len();
            Value::Cell(CellArray::from_dims(vec![l, 1], d))
        }
        other => other.clone(),
    };
    let gi = it.call_function("grp2idx", vec![gsel], 1)?;
    let Value::Num(gim) = &gi[0] else { return rt("grp2idx") };
    let ig: Vec<usize> = gim.data.iter().map(|v| *v as usize - 1).collect();
    let k = ig.iter().copied().max().map(|v| v + 1).unwrap_or(0);
    let mut cnt = vec![0.0; k];
    let mut sum = vec![0.0; k];
    for (v, &g) in x.iter().zip(&ig) {
        cnt[g] += 1.0;
        sum[g] += v;
    }
    let gmean: Vec<f64> = (0..k).map(|g| sum[g] / cnt[g]).collect();
    let mut ss = vec![0.0; k];
    for (v, &g) in x.iter().zip(&ig) {
        ss[g] += (v - gmean[g]).powi(2);
    }
    // grpstats reports the standard error of the mean; the variance is rebuilt from it.
    let gvar: Vec<f64> = (0..k)
        .map(|g| {
            let sd = if cnt[g] > 1.0 { (ss[g] / (cnt[g] - 1.0)).sqrt() } else { 0.0 };
            let sem = sd / cnt[g].sqrt();
            cnt[g] * sem * sem
        })
        .collect();
    let df: Vec<f64> = cnt.iter().map(|c| c - 1.0).collect();
    let sumdf: f64 = df.iter().sum();
    let vp = if sumdf > 0.0 { df.iter().zip(&gvar).map(|(d, v)| d * v).sum::<f64>() / sumdf } else { f64::NAN };
    let (p, bstat, bdf, statname): (f64, Value, Value, &str);
    if matches!(dorobust, 1 | 3 | 4 | 5 | 6) {
        let inc: Vec<usize> = (0..x.len()).filter(|&i| cnt[ig[i]] >= 2.0).collect();
        let nsp = cnt.iter().filter(|&&c| c < 2.0).count();
        let ngroups = k - nsp;
        let xt: Vec<f64> = match dorobust {
            1 | 6 => inc.iter().map(|&i| (x[i] - gmean[ig[i]]).powi(2)).collect(),
            3 => inc.iter().map(|&i| (x[i] - gmean[ig[i]]).abs()).collect(),
            4 => inc
                .iter()
                .map(|&i| {
                    let (c, gvr) = (cnt[ig[i]], gvar[ig[i]]);
                    let w = 0.5;
                    ((w + c - 2.0) * c * (x[i] - gmean[ig[i]]).powi(2) - w * (c - 1.0) * gvr) / ((c - 1.0) * (c - 2.0))
                })
                .collect(),
            _ => {
                let med: Vec<f64> = (0..k)
                    .map(|g| {
                        let mut v: Vec<f64> = x.iter().zip(&ig).filter(|(_, &gg)| gg == g).map(|(v, _)| *v).collect();
                        median_of(&mut v)
                    })
                    .collect();
                inc.iter().map(|&i| (x[i] - med[ig[i]]).abs()).collect()
            }
        };
        if ngroups > 1 {
            let gg: Vec<usize> = inc.iter().map(|&i| ig[i]).collect();
            let (f, d1, d2, pp) = oneway(it, &xt, &gg)?;
            p = pp;
            bstat = f;
            bdf = row(vec![d1, d2]);
        } else {
            p = f64::NAN;
            bstat = Value::num(f64::NAN);
            bdf = row(vec![0.0, xt.len() as f64 - ngroups as f64]);
        }
        statname = "fstat";
    } else {
        let pos: Vec<usize> = (0..k).filter(|&g| df[g] > 0.0).collect();
        let bd = (pos.len() as f64 - 1.0).max(0.0);
        let b = if bd > 0.0 && sumdf > 0.0 {
            let b0 = vp.ln() * sumdf - pos.iter().map(|&g| df[g] * gvar[g].ln()).sum::<f64>();
            let c = 1.0 + (pos.iter().map(|&g| 1.0 / df[g]).sum::<f64>() - 1.0 / sumdf) / (3.0 * bd);
            b0 / c
        } else {
            f64::NAN
        };
        let r = it.call_function("chi2cdf", vec![Value::num(b), Value::num(bd), Value::chars("upper")], 1)?;
        p = r[0].scalar()?;
        bstat = Value::num(b);
        bdf = Value::num(bd);
        statname = "chisqstat";
    }
    let single = is_single(xm);
    let to_single = |v: Value| match v {
        Value::Num(m) => Value::Num(finish(m, single)),
        other => other,
    };
    let mut out = vec![to_single(Value::num(p))];
    if nargout > 1 {
        let mut st = StructArray::scalar();
        st.set(statname, to_single(bstat));
        st.set("df", bdf);
        out.push(Value::Struct(st));
    }
    Ok(Some(out))
}

// ---------------------------------------------------------------------------------------------
// ansaribradley

/// Probabilities (below, equal, above) of the permutation distribution of the weighted count of
/// group 1 over the distinct values, against the observed statistic `tstar`.
fn ab_exact(n1: &[f64], n2: &[f64], wts: &[f64], tstar: f64) -> [f64; 3] {
    let m: usize = n1.iter().sum::<f64>() as usize;
    let seps = f64::EPSILON.sqrt();
    let mut states: std::collections::HashMap<(usize, i64), (f64, f64)> = std::collections::HashMap::new();
    states.insert((0, 0), (1.0, 0.0));
    for j in 0..wts.len() {
        let nj = n1[j] + n2[j];
        let mut next: std::collections::HashMap<(usize, i64), (f64, f64)> = std::collections::HashMap::new();
        for (&(c, _), &(pr, st)) in states.iter() {
            let maxx = (nj as usize).min(m - c);
            for xj in 0..=maxx {
                let pj = (crate::builtins_specfun::gammaln(nj + 1.0) - crate::builtins_specfun::gammaln(xj as f64 + 1.0) - crate::builtins_specfun::gammaln(nj - xj as f64 + 1.0)).exp();
                let ns = st + wts[j] * xj as f64;
                let key = ((c + xj), (ns / 1e-9).round() as i64);
                let e = next.entry(key).or_insert((0.0, ns));
                e.0 += pr * pj;
            }
        }
        states = next;
    }
    let mut p = [0.0; 3];
    let mut total = 0.0;
    for (&(c, _), &(pr, st)) in states.iter() {
        if c != m {
            continue;
        }
        total += pr;
        if st < tstar - seps {
            p[0] += pr;
        } else if st > tstar + seps {
            p[2] += pr;
        } else {
            p[1] += pr;
        }
    }
    [p[0] / total, p[1] / total, p[2] / total]
}

fn ansaribradley(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    if a.len() < 2 {
        return rt("Requires at least two input arguments");
    }
    let (Value::Num(xm), Value::Num(ym)) = (&a[0], &a[1]) else { return rt("Not enough input arguments.") };
    let mut alpha_v: Option<Value> = None;
    let mut tail_v: Option<Value> = None;
    let mut method_v: Option<Value> = None;
    let mut dim_v: Option<Value> = None;
    if a.len() >= 3 {
        if matches!(a[2], Value::Num(_)) && !a[2].is_text() {
            alpha_v = Some(a[2].clone());
            if a.len() >= 4 {
                tail_v = Some(a[3].clone());
                if a.len() >= 5 {
                    method_v = Some(a[4].clone());
                    if a.len() >= 6 {
                        dim_v = Some(a[5].clone());
                    }
                }
            }
        } else if a.len() == 3 {
            return rt("Alpha must be a scalar between 0 and 1.");
        } else {
            let rest = &a[2..];
            if rest.len() % 2 != 0 {
                return rt("Wrong number of arguments.");
            }
            for pair in rest.chunks(2) {
                let key = text_of(&pair[0]).to_lowercase();
                let names = ["alpha", "tail", "method", "dim"];
                let hits: Vec<usize> = names.iter().enumerate().filter(|(_, n)| !key.is_empty() && n.starts_with(&key)).map(|(i, _)| i).collect();
                if hits.len() != 1 {
                    return rt(format!("Invalid parameter name: {}.", text_of(&pair[0])));
                }
                let v = Some(pair[1].clone());
                match hits[0] {
                    0 => alpha_v = v,
                    1 => tail_v = v,
                    2 => method_v = v,
                    _ => dim_v = v,
                }
            }
        }
    }
    let alpha = match &alpha_v {
        None => 0.05,
        Some(v) if is_empty_val(v) => 0.05,
        Some(Value::Num(m)) if m.data.len() == 1 && m.data[0] > 0.0 && m.data[0] < 1.0 => m.data[0],
        _ => return rt("Alpha must be a scalar between 0 and 1."),
    };
    let bad_tail = || rt::<()>("Tail must be one of the character vectors or string scalars 'both', 'right', or 'left'.");
    let tail: f64 = match &tail_v {
        None => 0.0,
        Some(v) if is_empty_val(v) => 0.0,
        Some(v) if v.is_text() => {
            let t = v.text().unwrap_or_default().to_lowercase();
            let names = ["left", "both", "right"];
            match names.iter().position(|n| n.starts_with(&t)) {
                Some(i) if !t.is_empty() => i as f64 - 1.0,
                _ => return bad_tail().map(|_| None),
            }
        }
        Some(Value::Num(m)) if m.data.len() == 1 => m.data[0],
        _ => return bad_tail().map(|_| None),
    };
    let doexact: Option<bool> = match &method_v {
        None => None,
        Some(v) if is_empty_val(v) => None,
        Some(v) => {
            let hits: Vec<usize> = if v.is_text() {
                let t = v.text().unwrap_or_default().to_lowercase();
                ["on", "exact", "off", "approximate"].iter().enumerate().filter(|(_, n)| n.starts_with(&t)).map(|(i, _)| i + 1).collect()
            } else {
                vec![]
            };
            if hits.len() != 1 {
                return rt("The 'method' value must be either 'exact' or 'approximate'.");
            }
            Some(hits[0] == 1 || hits[0] == 2)
        }
    };
    let x_is_vec = is_vector_dims(&xm.dims);
    let y_is_vec = is_vector_dims(&ym.dims);
    let dim: usize = match &dim_v {
        Some(v) if !is_empty_val(v) => match v {
            Value::Num(m) if m.data.len() == 1 && (m.data[0] == 1.0 || m.data[0] == 2.0) => m.data[0] as usize,
            _ => return rt("Dimension must be an integer between 1 and 2."),
        },
        _ => xm.dims.iter().position(|&d| d != 1).map(|p| p + 1).unwrap_or(1),
    };
    // Columns of the samples along `dim`.
    let cols = |m: &Matrix, d: usize| -> Vec<Vec<f64>> {
        let (r, c) = (m.dims[0], m.dims[1]);
        if d == 1 {
            (0..c).map(|j| m.data[j * r..(j + 1) * r].to_vec()).collect()
        } else {
            (0..r).map(|i| (0..c).map(|j| m.data[j * r + i]).collect()).collect()
        }
    };
    let (xc, yc);
    if dim_v.as_ref().map(|v| is_empty_val(v)).unwrap_or(true) && x_is_vec && y_is_vec {
        // The vector y follows the orientation of x.
        xc = vec![xm.data.clone()];
        yc = vec![ym.data.clone()];
    } else {
        xc = cols(xm, dim);
        yc = cols(ym, dim);
        let (sx, sy) = (xm.dims.clone(), ym.dims.clone());
        let other = |d: &Vec<usize>| d[2 - dim];
        if other(&sx) != other(&sy) {
            return rt("The sizes of X and Y must be commensurate.");
        }
    }
    let nsamples = xc.len();
    let single = is_single(xm) || is_single(ym);
    let (mut ws, mut wstars, mut ps) = (vec![], vec![], vec![]);
    for j in 0..nsamples {
        let xv: Vec<f64> = xc[j].iter().copied().filter(|v| !v.is_nan()).collect();
        let yv: Vec<f64> = yc[j].iter().copied().filter(|v| !v.is_nan()).collect();
        let mut zs: Vec<(f64, u8)> = xv.iter().map(|&v| (v, 1u8)).chain(yv.iter().map(|&v| (v, 2u8))).collect();
        zs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let nn = zs.len();
        let sorted = Value::Num(Matrix::from_dims(vec![nn, 1], zs.iter().map(|z| z.0).collect()));
        let r = if nn == 0 {
            vec![]
        } else {
            let rr = it.call_function("tiedrank", vec![sorted, Value::num(0.0), Value::num(1.0)], 1)?;
            let Value::Num(rm) = &rr[0] else { return rt("tiedrank") };
            rm.data.clone()
        };
        let w: f64 = zs.iter().zip(&r).filter(|(z, _)| z.1 == 1).map(|(_, r)| *r).sum();
        let m = xv.len() as f64;
        let n = nn as f64 - m;
        let nf = nn as f64;
        let sumsq: f64 = r.iter().map(|v| v * v).sum();
        let (mean_w, std_w) = if nn % 2 == 0 {
            (m * (nf + 2.0) / 4.0, (m * n * (16.0 * sumsq - nf * (nf + 2.0).powi(2)) / (16.0 * nf * (nf - 1.0))).sqrt())
        } else {
            (m * (nf + 1.0).powi(2) / (4.0 * nf), (m * n * (16.0 * nf * sumsq - (nf + 1.0).powi(4)) / (16.0 * nf * nf * (nf - 1.0))).sqrt())
        };
        let wstar = if std_w > 0.0 {
            (w - mean_w) / std_w
        } else if w == mean_w {
            f64::NAN
        } else {
            (w - mean_w).signum() * f64::INFINITY
        };
        let exact = doexact.unwrap_or(xv.len() + yv.len() <= 25);
        let pvals: [f64; 3] = if m == 0.0 || n == 0.0 {
            [f64::NAN; 3]
        } else if exact {
            let mut u: Vec<f64> = zs.iter().map(|z| z.0).collect();
            u.dedup();
            let mut n1 = vec![0.0; u.len()];
            let mut n2 = vec![0.0; u.len()];
            let mut wts = vec![0.0; u.len()];
            for (z, rk) in zs.iter().zip(&r) {
                let k = u.partition_point(|&v| v < z.0);
                if z.1 == 1 {
                    n1[k] += 1.0;
                } else {
                    n2[k] += 1.0;
                }
                wts[k] = *rk;
            }
            ab_exact(&n1, &n2, &wts, w)
        } else {
            let pp = it.call_function("normcdf", vec![Value::num(-wstar.abs())], 1)?[0].scalar()?;
            if wstar < 0.0 { [pp, 0.0, 1.0 - pp] } else { [1.0 - pp, 0.0, pp] }
        };
        let p = if tail == 0.0 {
            (2.0 * (pvals[1] + pvals[0].min(pvals[2]))).min(1.0)
        } else if tail == 1.0 {
            pvals[1] + pvals[0]
        } else if tail == -1.0 {
            pvals[1] + pvals[2]
        } else {
            0.0
        };
        ws.push(w);
        wstars.push(wstar);
        ps.push(p);
    }
    let shape = if dim == 1 { vec![1, nsamples] } else { vec![nsamples, 1] };
    let mk = |v: Vec<f64>| Value::Num(finish(Matrix::from_dims(shape.clone(), v), single));
    let h: Vec<f64> = ps.iter().map(|&p| if p.is_nan() { f64::NAN } else if p <= alpha { 1.0 } else { 0.0 }).collect();
    let mut out = vec![mk(h), mk(ps)];
    if nargout > 2 {
        let mut st = StructArray::scalar();
        st.set("W", mk(ws));
        st.set("Wstar", mk(wstars));
        out.push(Value::Struct(st));
    }
    Ok(Some(out))
}

// ---------------------------------------------------------------------------------------------
// lillietest

const LILLIE_ALPHAS: [f64; 12] = [0.001, 0.0015, 0.002, 0.005, 0.01, 0.015, 0.02, 0.05, 0.10, 0.15, 0.20, 0.50];

fn lilliefors_cvs(distr: usize, n: usize) -> Vec<f64> {
    let (q, c) = match distr {
        0 => (&LILLIE_Q1, &LILLIE_C1),
        1 => (&LILLIE_Q2, &LILLIE_C2),
        _ => (&LILLIE_Q3, &LILLIE_C3),
    };
    if n <= 20 {
        q[n - 4].to_vec()
    } else {
        let nf = n as f64;
        (0..12).map(|j| c[0][j] / nf.sqrt() - c[1][j] / nf - c[2][j] / nf.powf(1.5)).collect()
    }
}

/// The Kolmogorov-Smirnov distance between the sample and the distribution fitted to it by
/// maximum likelihood (0 normal, 1 exponential, 2 extreme value).
fn lillie_stat(it: &mut Interpreter, x: &[f64], distr: usize) -> Res<f64> {
    let n = x.len() as f64;
    let mut s = x.to_vec();
    s.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mut u: Vec<f64> = Vec::new();
    let mut cum: Vec<f64> = Vec::new();
    for (i, v) in s.iter().enumerate() {
        if u.last() == Some(v) {
            *cum.last_mut().unwrap() = (i + 1) as f64 / n;
        } else {
            u.push(*v);
            cum.push((i + 1) as f64 / n);
        }
    }
    let mean = x.iter().sum::<f64>() / n;
    let uv = Value::Num(Matrix::from_dims(vec![u.len(), 1], u.clone()));
    let null = match distr {
        0 => {
            let sd = (x.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / (n - 1.0)).sqrt();
            it.call_function("normcdf", vec![uv, Value::num(mean), Value::num(sd)], 1)?
        }
        1 => it.call_function("expcdf", vec![uv, Value::num(mean)], 1)?,
        _ => {
            let ph = it.call_function("evfit", vec![Value::Num(Matrix::from_dims(vec![x.len(), 1], x.to_vec()))], 1)?;
            let Value::Num(pm) = &ph[0] else { return rt("evfit") };
            it.call_function("evcdf", vec![uv, Value::num(pm.data[0]), Value::num(pm.data[1])], 1)?
        }
    };
    let Value::Num(nm) = &null[0] else { return rt("cdf") };
    if nm.data.iter().any(|v| v.is_nan()) {
        return Ok(f64::NAN);
    }
    let mut ks = f64::NEG_INFINITY;
    for i in 0..u.len() {
        let lo = if i == 0 { 0.0 } else { cum[i - 1] };
        ks = ks.max((lo - nm.data[i]).abs()).max((cum[i] - nm.data[i]).abs());
    }
    Ok(ks)
}

fn lillie_mc(it: &mut Interpreter, ks: f64, n: usize, alpha: f64, distr: usize, tol: f64) -> Res<(f64, f64)> {
    let vartol = tol * tol;
    let (mut crit, mut p, mut tot) = (0.0, 0.0, 0.0f64);
    let mut reps_min = 1000.0f64;
    loop {
        let old = tot;
        let reps = (reps_min - old).ceil().max(1.0) as usize;
        let mut stats = Vec::with_capacity(reps);
        for _ in 0..reps {
            let gen = if distr == 0 { "randn" } else { "rand" };
            let r = it.call_function(gen, vec![Value::num(n as f64), Value::num(1.0)], 1)?;
            let Value::Num(m) = &r[0] else { return rt("rand") };
            let sample: Vec<f64> = match distr {
                0 => m.data.clone(),
                1 => m.data.iter().map(|u| -(1.0 - u).ln()).collect(),
                _ => m.data.iter().map(|u| (-(1.0 - u).ln()).ln()).collect(),
            };
            stats.push(lillie_stat(it, &sample, distr)?);
        }
        let q = it.call_function("prctile", vec![Value::Num(Matrix::from_dims(vec![reps, 1], stats.clone())), Value::num(100.0 * (1.0 - alpha))], 1)?;
        let crit_mc = q[0].scalar()?;
        let p_mc = stats.iter().filter(|&&s| s > ks).count() as f64 / reps as f64;
        tot = old + reps as f64;
        crit = (old * crit + reps as f64 * crit_mc) / tot;
        p = (old * p + reps as f64 * p_mc) / tot;
        let sepsq = (p * (1.0 - p) / tot).max(1.0 / (tot * tot));
        if sepsq < vartol {
            break;
        }
        reps_min = 1.2 * (tot * sepsq) / vartol;
    }
    if ks.is_nan() {
        p = 0.0;
    }
    Ok((crit, p))
}

fn lillietest(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    if a.is_empty() {
        return rt("Not enough input arguments.");
    }
    let Some(Value::Num(xm)) = a.first() else { return rt("Input sample X must be a vector of real values.") };
    if !is_vector_dims(&xm.dims) || xm.is_complex() {
        return rt("Input sample X must be a vector of real values.");
    }
    let x: Vec<f64> = xm.data.iter().copied().filter(|v| !v.is_nan()).collect();
    if x.len() < 4 {
        return rt("Sample vector X must have at least 4 valid observations.");
    }
    let mut alpha_v: Option<Value> = None;
    let mut distr_v: Option<Value> = None;
    let mut mctol_v: Option<Value> = None;
    if a.len() > 1 {
        let first = &a[1];
        if (matches!(first, Value::Num(_)) && !first.is_text()) || is_empty_val(first) {
            alpha_v = Some(first.clone());
            if a.len() >= 3 {
                distr_v = Some(a[2].clone());
                if a.len() == 4 {
                    mctol_v = Some(a[3].clone());
                }
            }
        } else {
            if (a.len() - 1) % 2 != 0 {
                return rt("Wrong number of arguments.");
            }
            for pair in a[1..].chunks(2) {
                let key = text_of(&pair[0]).to_lowercase();
                let names = ["distribution", "alpha", "mctol"];
                let hits: Vec<usize> = names.iter().enumerate().filter(|(_, n)| !key.is_empty() && n.starts_with(&key)).map(|(i, _)| i).collect();
                if hits.len() != 1 {
                    return rt(format!("Invalid parameter name: {}.", text_of(&pair[0])));
                }
                let v = Some(pair[1].clone());
                match hits[0] {
                    0 => distr_v = v,
                    1 => alpha_v = v,
                    _ => mctol_v = v,
                }
            }
        }
    }
    let alpha = match &alpha_v {
        None => 0.05,
        Some(v) if is_empty_val(v) => 0.05,
        Some(Value::Num(m)) if m.data.len() == 1 && 0.0 < m.data[0] && m.data[0] < 1.0 => m.data[0],
        _ => return rt("Significance level ALPHA must be a scalar between 0 and 1."),
    };
    let names = ["normal", "exponential", "ev", "extreme value"];
    let distr: usize = match &distr_v {
        None => 0,
        Some(v) if is_empty_val(v) => 0,
        Some(v) if v.is_text() => {
            let k = param_val(v, &names, "distribution name", true)?;
            if k == 3 { 2 } else { k }
        }
        _ => return rt("DISTR must be 'norm', 'exp', or 'ev'."),
    };
    let mctol = match &mctol_v {
        None => None,
        Some(v) if is_empty_val(v) => None,
        Some(Value::Num(m)) if m.data.len() == 1 && m.data[0] > 0.0 => Some(m.data[0]),
        _ => return rt("Monte-Carlo standard error tolerance MCTOL must be a positive scalar value."),
    };
    let n = x.len();
    let ks = lillie_stat(it, &x, distr)?;
    let (cv, p);
    let mut from_root = false;
    if let Some(tol) = mctol {
        let (c, pp) = lillie_mc(it, ks, n, alpha, distr, tol)?;
        cv = c;
        p = pp;
    } else {
        if alpha < LILLIE_ALPHAS[0] || LILLIE_ALPHAS[11] < alpha {
            return rt("Significance level ALPHA is outside the range of the tabulated values.\nUse a value in the interval [0.001, 0.5], or use the MCTOL input argument.");
        }
        let cvs = lilliefors_cvs(distr, n);
        cv = pchip_at(&LILLIE_ALPHAS, &cvs, alpha);
        let mut pv = f64::NAN;
        if nargout > 1 {
            if ks < cvs[11] {
                it.emit_warning("P is greater than the largest tabulated value, returning 0.5.".into(), "stats:lillietest:OutOfRangePHigh".into());
                pv = LILLIE_ALPHAS[11];
            } else if cvs[0] < ks {
                it.emit_warning("P is less than the smallest tabulated value, returning 0.001.".into(), "stats:lillietest:OutOfRangePLow".into());
                pv = LILLIE_ALPHAS[0];
            } else if ks == cvs[0] {
                pv = LILLIE_ALPHAS[0];
            } else if ks.is_nan() {
                pv = 0.0;
            } else {
                let i = cvs.iter().position(|&c| ks > c).unwrap();
                let (mut lo, mut hi) = (LILLIE_ALPHAS[i - 1], LILLIE_ALPHAS[i]);
                let f = |t: f64| pchip_at(&LILLIE_ALPHAS, &cvs, t) - ks;
                let flo = f(lo);
                for _ in 0..200 {
                    let mid = 0.5 * (lo + hi);
                    if (f(mid) > 0.0) == (flo > 0.0) {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                pv = 0.5 * (lo + hi);
                from_root = true;
            }
        }
        p = pv;
    }
    let h = if ks > cv { 1.0 } else { 0.0 } + if ks.is_nan() { 1.0 } else { 0.0 };
    let single = is_single(xm);
    let sv = |v: f64| Value::Num(finish(Matrix::scalar(v), single));
    Ok(Some(vec![Value::num(h), if from_root { sv(p) } else { Value::num(p) }, sv(ks), Value::num(cv)]))
}
