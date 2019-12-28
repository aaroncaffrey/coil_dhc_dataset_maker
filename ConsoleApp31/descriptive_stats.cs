using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Runtime.InteropServices;

namespace Coil_DHC_Dataset_Maker
{

    public class descriptive_stats
    {
        public string member_id_name;
        public uint count;
        public uint count_zero_values;
        public uint count_non_zero_values;
        public uint count_distinct_values;
        public double sum;
        public double mean_arithmetic;
        //public double mean_geometric;
        //public double mean_harmonic;
        public double min;
        public double max;
        public double range;
        public double mid_range;
        public double variance;
        public double dev_standard;
        public double skewness;
        public double kurtosis;
        public double interquartile_range;
        public double median_q1;
        public double median_q2;
        public double median_q3;
        public double sum_of_error;
        public double sum_of_error_square;
        public double mode;
        public double mad_mean_arithmetic;
        //public double mad_mean_harmonic;
        //public double mad_mean_geometric;
        public double mad_median_q1;
        public double mad_median_q2;
        public double mad_median_q3;
        public double mad_mode;
        public double mad_mid_range;

        public descriptive_stats interquartile_range_descriptive_stats;
        public descriptive_stats intervals_descriptive_stats;

        public class descriptive_stats_encoding_options
        {
            public string options_name = "";

            public bool count = false;
            public bool count_zero_values = false;
            public bool count_non_zero_values = false;
            public bool count_distinct_values = false;
            public bool sum = false;
            public bool mean_arithmetic = false;
            public bool min = false;
            public bool max = false;
            public bool range = false;
            public bool mid_range = false;
            public bool variance = false;
            public bool dev_standard = false;
            public bool skewness = false;
            public bool kurtosis = false;
            public bool interquartile_range = false;
            public bool median_q1 = false;
            public bool median_q2 = false;
            public bool median_q3 = false;
            public bool sum_of_error = false;
            public bool sum_of_error_square = false;
            public bool mode = false;
            public bool mad_mean_arithmetic = false;
            public bool mad_median_q1 = false;
            public bool mad_median_q2 = false;
            public bool mad_median_q3 = false;
            public bool mad_mode = false;
            public bool mad_mid_range = false;

            public descriptive_stats_encoding_options(string name, bool enable = false)
            {
                options_name = name;

                count = enable;
                count_zero_values = enable;
                count_non_zero_values = enable;
                count_distinct_values = enable;
                sum = enable;
                mean_arithmetic = enable;
                min = enable;
                max = enable;
                range = enable;
                mid_range = enable;
                variance = enable;
                dev_standard = enable;
                skewness = enable;
                kurtosis = enable;
                interquartile_range = enable;
                median_q1 = enable;
                median_q2 = enable;
                median_q3 = enable;
                sum_of_error = enable;
                sum_of_error_square = enable;
                mode = enable;
                mad_mean_arithmetic = enable;
                mad_median_q1 = enable;
                mad_median_q2 = enable;
                mad_median_q3 = enable;
                mad_mode = enable;
                mad_mid_range = enable;
            }

            public static descriptive_stats_encoding_options options_average()
            {
                return new descriptive_stats_encoding_options("mean", false)
                {
                    mean_arithmetic = true,
                };
            }

            public static descriptive_stats_encoding_options options_average_sd()
            {
                return new descriptive_stats_encoding_options("mean", false)
                {
                    mean_arithmetic = true,
                    dev_standard = true,
                };
            }
        }

        public static readonly descriptive_stats_encoding_options options_default = descriptive_stats_encoding_options.options_average();

        public static readonly descriptive_stats_encoding_options options_all = new descriptive_stats_encoding_options("all", true);


        public List<(string member_id, string perspective_id, double perspective_value)> encode(descriptive_stats_encoding_options descriptive_stats_encoding_options = null, bool encode_intervals_descriptive_stats = false, bool encode_interquartile_range_descriptive_stats = false)
        {
            return descriptive_stats.encode(this, descriptive_stats_encoding_options, encode_intervals_descriptive_stats, encode_interquartile_range_descriptive_stats);
        }

        public static List<(string member_id, string perspective_id, double perspective_value)> encode(descriptive_stats stats, descriptive_stats_encoding_options descriptive_stats_encoding_options = null, bool encode_intervals_descriptive_stats = false, bool encode_interquartile_range_descriptive_stats = false)
        {
            var result = new List<(string member_id, string perspective_id, double perspective_value)>();

            if (stats == null) stats = new descriptive_stats(null, "");

            if (encode_intervals_descriptive_stats && stats.intervals_descriptive_stats != null)
            {
                var x = encode(stats.intervals_descriptive_stats, descriptive_stats_encoding_options, false, false);

                result.AddRange(x);
            }

            if (encode_interquartile_range_descriptive_stats && stats.interquartile_range_descriptive_stats != null)
            {
                var y = encode(stats.interquartile_range_descriptive_stats, descriptive_stats_encoding_options, false, false);

                result.AddRange(y);
            }

            if (descriptive_stats_encoding_options == null)
            {
                descriptive_stats_encoding_options = options_default;
            }

            var z = new List<(string member_id, string perspective_id, double perspective_value)>();
            if (descriptive_stats_encoding_options.count) { z.Add((member_id: $"{stats.member_id_name}", perspective_id: $"{nameof(stats.count)}", (double)stats.count)); }
            if (descriptive_stats_encoding_options.count_zero_values) { z.Add((member_id: $"{stats.member_id_name}", perspective_id: $"{nameof(stats.count_zero_values)}", (double)stats.count_zero_values)); }
            if (descriptive_stats_encoding_options.count_non_zero_values) { z.Add((member_id: $"{stats.member_id_name}", perspective_id: $"{nameof(stats.count_non_zero_values)}", (double)stats.count_non_zero_values)); }
            if (descriptive_stats_encoding_options.count_distinct_values) { z.Add((member_id: $"{stats.member_id_name}", perspective_id: $"{nameof(stats.count_distinct_values)}", (double)stats.count_distinct_values)); }
            if (descriptive_stats_encoding_options.min) { z.Add((member_id: $"{stats.member_id_name}", perspective_id: $"{nameof(stats.min)}", stats.min)); }
            if (descriptive_stats_encoding_options.max) { z.Add((member_id: $"{stats.member_id_name}", perspective_id: $"{nameof(stats.max)}", stats.max)); }
            if (descriptive_stats_encoding_options.range) { z.Add((member_id: $"{stats.member_id_name}", perspective_id: $"{nameof(stats.range)}", stats.range)); }
            if (descriptive_stats_encoding_options.sum) { z.Add((member_id: $"{stats.member_id_name}", perspective_id: $"{nameof(stats.sum)}", stats.sum)); }
            if (descriptive_stats_encoding_options.mid_range) { z.Add((member_id: $"{stats.member_id_name}", perspective_id: $"{nameof(stats.mid_range)}", stats.mid_range)); }
            if (descriptive_stats_encoding_options.mode) { z.Add((member_id: $"{stats.member_id_name}", perspective_id: $"{nameof(stats.mode)}", stats.mode)); }
            if (descriptive_stats_encoding_options.median_q1) { z.Add((member_id: $"{stats.member_id_name}", perspective_id: $"{nameof(stats.median_q1)}", stats.median_q1)); }
            if (descriptive_stats_encoding_options.median_q2) { z.Add((member_id: $"{stats.member_id_name}", perspective_id: $"{nameof(stats.median_q2)}", stats.median_q2)); }
            if (descriptive_stats_encoding_options.median_q3) { z.Add((member_id: $"{stats.member_id_name}", perspective_id: $"{nameof(stats.median_q3)}", stats.median_q3)); }
            if (descriptive_stats_encoding_options.mean_arithmetic) { z.Add((member_id: $"{stats.member_id_name}", perspective_id: $"{nameof(stats.mean_arithmetic)}", stats.mean_arithmetic)); }
            //if (descriptive_stats_encoding_options.mean_harmonic) { z.Add((member_id: $"{stats.name}", perspective_id: $"{nameof(stats.mean_harmonic)}",stats.mean_harmonic)); }
            //if (descriptive_stats_encoding_options.mean_geometric) { z.Add((member_id: $"{stats.name}", perspective_id: $"{nameof(stats.mean_geometric)}",stats.mean_geometric)); }
            if (descriptive_stats_encoding_options.variance) { z.Add((member_id: $"{stats.member_id_name}", perspective_id: $"{nameof(stats.variance)}", stats.variance)); }
            if (descriptive_stats_encoding_options.dev_standard) { z.Add((member_id: $"{stats.member_id_name}", perspective_id: $"{nameof(stats.dev_standard)}", stats.dev_standard)); }
            if (descriptive_stats_encoding_options.mad_mean_arithmetic) { z.Add((member_id: $"{stats.member_id_name}", perspective_id: $"{nameof(stats.mad_mean_arithmetic)}", stats.mad_mean_arithmetic)); }
            //if (descriptive_stats_encoding_options.mad_mean_harmonic) { z.Add((member_id: $"{stats.name}", perspective_id: $"{nameof(stats.mad_mean_harmonic)}",stats.mad_mean_harmonic)); }
            //if (descriptive_stats_encoding_options.mad_mean_geometric) { z.Add((member_id: $"{stats.name}", perspective_id: $"{nameof(stats.mad_mean_geometric)}",stats.mad_mean_geometric)); }
            if (descriptive_stats_encoding_options.mad_median_q1) { z.Add((member_id: $"{stats.member_id_name}", perspective_id: $"{nameof(stats.mad_median_q1)}", stats.mad_median_q1)); }
            if (descriptive_stats_encoding_options.mad_median_q2) { z.Add((member_id: $"{stats.member_id_name}", perspective_id: $"{nameof(stats.mad_median_q2)}", stats.mad_median_q2)); }
            if (descriptive_stats_encoding_options.mad_median_q3) { z.Add((member_id: $"{stats.member_id_name}", perspective_id: $"{nameof(stats.mad_median_q3)}", stats.mad_median_q3)); }
            if (descriptive_stats_encoding_options.mad_mode) { z.Add((member_id: $"{stats.member_id_name}", perspective_id: $"{nameof(stats.mad_mode)}", stats.mad_mode)); }
            if (descriptive_stats_encoding_options.mad_mid_range) { z.Add((member_id: $"{stats.member_id_name}", perspective_id: $"{nameof(stats.mad_mid_range)}", stats.mad_mid_range)); }
            if (descriptive_stats_encoding_options.sum_of_error) { z.Add((member_id: $"{stats.member_id_name}", perspective_id: $"{nameof(stats.sum_of_error)}", stats.sum_of_error)); }
            if (descriptive_stats_encoding_options.sum_of_error_square) { z.Add((member_id: $"{stats.member_id_name}", perspective_id: $"{nameof(stats.sum_of_error_square)}", stats.sum_of_error_square)); }
            if (descriptive_stats_encoding_options.interquartile_range) { z.Add((member_id: $"{stats.member_id_name}", perspective_id: $"{nameof(stats.interquartile_range)}", stats.interquartile_range)); }
            if (descriptive_stats_encoding_options.skewness) { z.Add((member_id: $"{stats.member_id_name}", perspective_id: $"{nameof(stats.skewness)}", stats.skewness)); }
            if (descriptive_stats_encoding_options.kurtosis) { z.Add((member_id: $"{stats.member_id_name}", perspective_id: $"{nameof(stats.kurtosis)}", stats.kurtosis)); }
            result.AddRange(z);


            return result;
        }

        public static descriptive_stats get_stat_values(double[] data, string member_id_name, bool preserve_values = false, bool presorted = false, bool calc_interquartile_range_stats = false, bool calc_interval_stats = false)
        {
            return new descriptive_stats(data, member_id_name, preserve_values, presorted, calc_interquartile_range_stats, calc_interval_stats);
        }

        public descriptive_stats(double[] data, string member_id_name, bool preserve_values = false, bool presorted = false, bool calc_interquartile_range_stats = false, bool calc_interval_stats = false)
        {
            this.member_id_name = member_id_name;

            if (data == null || data.Length == 0 || data.All(a => a == 0))
            {
                if (calc_interquartile_range_stats)
                {
                    this.interquartile_range_descriptive_stats = new descriptive_stats(new double[0], $"{this.member_id_name}_{nameof(this.interquartile_range_descriptive_stats)}", preserve_values: false, presorted: true, calc_interquartile_range_stats: false, calc_interval_stats: false);
                }

                if (calc_interval_stats)
                {
                    this.intervals_descriptive_stats = new descriptive_stats(new double[0], $"{this.member_id_name}_{nameof(this.intervals_descriptive_stats)}", preserve_values: false, presorted: true, calc_interquartile_range_stats: false, calc_interval_stats: false);
                }

                return;
            }

            double[] sorted_data;

            if (presorted)
            {
                sorted_data = data;
            }
            else
            {
                sorted_data = new double[data.Length];
                data.CopyTo(sorted_data, 0);
                Array.Sort(sorted_data);
            }

            this.count_distinct_values = (uint)sorted_data.Distinct().Count();

            //var cum_product = (double)1.0; // to calculate geometric mean
            //var cum_reciprocal = (double)0.0; // to calculate harmonic mean

            // First iteration
            for (var i = 0; i < sorted_data.Length; i++)
            {
                if (sorted_data[i] == (double)0.0)
                {
                    this.count_zero_values++;
                    continue;
                }
                this.count_non_zero_values++;

                //cum_product *= sorted_data[i];
                //cum_reciprocal += ((double)1.0) / sorted_data[i];
            }

            this.count = (uint)sorted_data.Length;
            var n = (double)this.count; // use a shorter variable in double type
            this.sum = sorted_data.Sum();
            this.mean_arithmetic = this.sum / n;
            //this.mean_geometric = (double)Math.Pow((double)cum_product, (double)((double)1.0 / n));
            //this.mean_harmonic = ((double)1.0) / (cum_reciprocal / n); // see http://mathworld.wolfram.com/HarmonicMean.html
            //this.range = this.max - this.min;

            // second loop, calculate Stdev, sum of errors
            //double[] eSquares = new double[data.Length];
            var m1 = (double)0.0;
            var m2 = (double)0.0;
            var m3 = (double)0.0; // for skewness calculation
            var m4 = (double)0.0; // for kurtosis calculation
                                  // for skewness
            for (int i = 0; i < sorted_data.Length; i++)
            {
                var m = sorted_data[i] - this.mean_arithmetic;
                var mPow2 = m * m;
                var mPow3 = mPow2 * m;
                var mPow4 = mPow3 * m;

                m1 += Math.Abs(m);

                m2 += mPow2;

                // calculate skewness
                m3 += mPow3;

                // calculate skewness
                m4 += mPow4;

            }

            this.sum_of_error = m1;
            this.sum_of_error_square = m2;

            this.variance = this.count <= 1 ? (double)0.0 : this.sum_of_error_square / ((double)this.count - 1);
            this.dev_standard = this.variance == (double)0.0 ? (double)0.0 : (double)Math.Sqrt((double)this.variance);

            var skew_cum = (double)0.0;
            if (this.dev_standard != (double)0.0)
            {
                for (int i = 0; i < sorted_data.Length; i++)
                {
                    skew_cum += (double)Math.Pow((double)((sorted_data[i] - this.mean_arithmetic) / this.dev_standard), 3);
                }
            }

            if (skew_cum != 0)
            {
                var x = (n - 1);
                var y = (n - 2);
                if (x != 0 && y != 0)// && skew_cum != 0)
                {
                    this.skewness = n / x / y * skew_cum;
                }
            }


            // kurtosis: see http://en.wikipedia.org/wiki/Kurtosis (heading: Sample Kurtosis)
            var m2_2 = (double)Math.Pow((double)this.sum_of_error_square, 2);
            if (m2_2 != (double)0.0)
            {
                var x = ((n - 2) * (n - 3));
                var y = (m4 / m2_2);
                var z = ((double)Math.Pow((double)(n - 1), 2));

                if (x != 0 && y != 0 && z != 0)
                {
                    this.kurtosis = ((n + 1) * n * (n - 1)) / x * y - 3 * z / x;
                }
            }
            // calculate quartiles


            this.min = sorted_data[0];
            this.max = sorted_data[sorted_data.Length - 1];
            this.range = this.max - this.min;
            this.mid_range = (this.max + this.min) / (double)2.0;

            this.median_q1 = percentile(sorted_data, 25);
            this.median_q3 = percentile(sorted_data, 75);
            this.median_q2 = percentile(sorted_data, 50);
            this.interquartile_range = Math.Abs(this.median_q3 - this.median_q1);

            var sorted_data_groups = sorted_data.GroupBy(x => x).ToArray();
            var sorted_data_max_group_count = sorted_data_groups.Max(g => g.Count());
            var modes = sorted_data_groups.Where(a => a.Count() == sorted_data_max_group_count).ToList();
            var mode_mean = modes.Select(a => a.Key).DefaultIfEmpty(0).Average();
            this.mode = mode_mean;

            this.mad_mean_arithmetic = mad(sorted_data, this.mean_arithmetic);
            //this.mad_mean_geometric = mad(sorted_data, this.mean_geometric);
            //this.mad_mean_harmonic = mad(sorted_data, this.mean_harmonic);
            this.mad_median_q1 = mad(sorted_data, this.median_q1);
            this.mad_median_q2 = mad(sorted_data, this.median_q2);
            this.mad_median_q3 = mad(sorted_data, this.median_q3);
            this.mad_mid_range = mad(sorted_data, this.mid_range);
            this.mad_mode = mad(sorted_data, this.mode);



            fix_double(nameof(sum), ref sum);
            fix_double(nameof(mean_arithmetic), ref mean_arithmetic);
            //fix_double(nameof(mean_geometric), ref mean_geometric);
            //fix_double(nameof(mean_harmonic), ref mean_harmonic);
            fix_double(nameof(min), ref min);
            fix_double(nameof(max), ref max);
            fix_double(nameof(range), ref range);
            fix_double(nameof(mid_range), ref mid_range);
            fix_double(nameof(variance), ref variance);
            fix_double(nameof(dev_standard), ref dev_standard);
            fix_double(nameof(skewness), ref skewness);
            fix_double(nameof(kurtosis), ref kurtosis);
            fix_double(nameof(interquartile_range), ref interquartile_range);
            fix_double(nameof(median_q1), ref median_q1);
            fix_double(nameof(median_q2), ref median_q2);
            fix_double(nameof(median_q3), ref median_q3);
            fix_double(nameof(sum_of_error), ref sum_of_error);
            fix_double(nameof(sum_of_error_square), ref sum_of_error_square);
            fix_double(nameof(mode), ref mode);
            fix_double(nameof(mad_mean_arithmetic), ref mad_mean_arithmetic);
            //fix_double(nameof(mad_mean_harmonic), ref mad_mean_harmonic);
            //fix_double(nameof(mad_mean_geometric), ref mad_mean_geometric);
            fix_double(nameof(mad_median_q1), ref mad_median_q1);
            fix_double(nameof(mad_median_q2), ref mad_median_q2);
            fix_double(nameof(mad_median_q3), ref mad_median_q3);
            fix_double(nameof(mad_mode), ref mad_mode);
            fix_double(nameof(mad_mid_range), ref mad_mid_range);


            if (calc_interquartile_range_stats)
            {
                var interquartile_data = sorted_data.Where(a => a >= this.median_q1 && a <= this.median_q3).ToArray();


                this.interquartile_range_descriptive_stats = new descriptive_stats(interquartile_data, $"{this.member_id_name}_{nameof(this.interquartile_range_descriptive_stats)}", preserve_values: false, presorted: true,
                    calc_interquartile_range_stats: false, calc_interval_stats: false);
            }

            if (calc_interval_stats)
            {
                var data_point_distances = new double[sorted_data.Length - 1];

                for (var i = 1; i < sorted_data.Length; i++)
                {
                    data_point_distances[i - 1] = Math.Abs(sorted_data[i] - sorted_data[i - 1]);
                }
                Array.Sort(data_point_distances);

                this.intervals_descriptive_stats = new descriptive_stats(data_point_distances, $"{this.member_id_name}_{nameof(this.intervals_descriptive_stats)}", preserve_values: false, presorted: true,
                    calc_interquartile_range_stats: false, calc_interval_stats: false);
            }

            if (!preserve_values)
            {
                data = null;
                sorted_data = null;
            }
        }

        public static void fix_double(string name, ref double value)
        {
            const double c_double_max = (double)1.79769e+308;
            const double c_double_min = (double)-c_double_max;
            const double double_zero = (double)0;

            var output = false;

            if (double.IsPositiveInfinity(value) || value >= c_double_max || value >= double.MaxValue)
            {
                if (output) Console.WriteLine("fix_double: " + name + " = " + value.ToString("G17") + " is positive infinity.");

                value = c_double_max;
            }
            else if (double.IsNegativeInfinity(value) || value <= c_double_min || value <= double.MinValue)
            {
                if (output) Console.WriteLine("fix_double: " + name + " = " + value.ToString("G17") + " is negative infinity.");

                value = c_double_min;
            }
            else if (double.IsNaN(value))
            {
                if (output) Console.WriteLine("fix_double: " + name + " = " + value.ToString("G17") + " is not a number.");

                value = double_zero;
            }
        }

        public static double fix_double(double value)
        {
            const double c_double_max = (double)1.79769e+308;
            const double c_double_min = (double)-c_double_max;
            const double double_zero = (double)0;

            if (double.IsPositiveInfinity(value) || value >= c_double_max || value >= double.MaxValue)
            {
                value = c_double_max;
            }
            else if (double.IsNegativeInfinity(value) || value <= c_double_min || value <= double.MinValue)
            {
                value = c_double_min;
            }
            else if (double.IsNaN(value))
            {
                value = double_zero;
            }

            return value;
        }

        public static double mad(double[] values, double? centre)
        {
            if (values.Length == 0) return 0;

            if (centre == null) centre = values.Average();

            var mad = values.Sum(a => Math.Abs(centre.Value - a)) / (double)values.Length;

            return mad;
        }

        public static double percentile(double[] sorted_data, double p)
        {
            if (sorted_data.Length == 0) return 0;
            if (sorted_data.Length == 1) return sorted_data[0];

            // algo derived from Aczel pg 15 bottom
            if (p >= (double)100.0)
            {
                return sorted_data[sorted_data.Length - 1];
            }

            var position = (double)(sorted_data.Length + 1) * p / ((double)100.0);
            var left_number = (double)0.0;
            var right_number = (double)0.0;

            double n = p / ((double)100.0) * (sorted_data.Length - 1) + ((double)1.0);

            if (position >= 1)
            {
                left_number = sorted_data[(int)System.Math.Floor(n) - 1];
                right_number = sorted_data[(int)System.Math.Floor(n)];
            }
            else
            {
                left_number = sorted_data[0];
                right_number = sorted_data[1];
            }

            if (left_number == right_number)
            {
                return left_number;
            }
            else
            {
                double part = n - Math.Floor(n);
                return left_number + part * (right_number - left_number);
            }
        }



    }
}
