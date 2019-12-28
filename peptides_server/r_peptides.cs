using System;
using System.Collections.Generic;
using System.Linq;
using System.Security.Cryptography.X509Certificates;
using System.Text;
using System.Threading.Tasks;
using RDotNet;

namespace Coil_DHC_Dataset_Maker
{
    public static class r_peptides
    {
        public static readonly object engine_lock = new object();

        //public static REngine engine = null;//init_r();

        private static uint _key = 1;
        private static readonly object _key_lock = new object();

        public static int id;

        public static string Key
        {
            get
            {
                lock (_key_lock)
                {
                    _key++;
                    return $"{id}_{_key}";
                }
            }
        }

        //public r_peptides()
        //{
        //    engine = init_r();
        //}

        private static bool need_init = true;

        public static REngine init_r()
        {

            StartupParameter rinit = new StartupParameter();

            rinit.Quiet = true;
            rinit.Interactive = false;

            rinit.RHome = $@"C:\Program Files\R\R-3.4.4\";

            //REngine.SetEnvironmentVariables();
            

            var engine1 = REngine.GetInstance(@"C:\Program Files\R\R-3.4.4\bin\x64\R.dll", true, rinit);

            if (need_init)
            {
                need_init = false;
                var r_init_cmds = $@"
            #install.packages(""devtools"", repos=""http://cran.us.r-project.org"")
            #library(devtools)
            #install_github(""dosorio/Peptides"")
            sink(""NUL"")
            library(Peptides)
            if (!exists(""AAdata"")) data(AAdata)
            ";
                r_init_cmds.Split(new char[] {'\r', '\n'}).Where(a => !string.IsNullOrWhiteSpace(a) && !a.Trim().StartsWith("#")).ToList().ForEach(a => engine1.Evaluate(a));

            }

            return engine1;
        }

        public static List<subsequence_classification_data.feature_info> calculate_r_peptides_classification_data_template = null;


        public static object get_values_lock = new object();

        public static List<subsequence_classification_data.feature_info> get_values(string seq, int id) //List<(string name, double value)>
        {
            lock (get_values_lock)
            {
                r_peptides.id = id;

#if DEBUG
                if (Program.verbose_debug) Program.WriteLine($"{nameof(get_values)}(string seq);");
#endif

                var engine = init_r();
                {



                    if (string.IsNullOrWhiteSpace(seq))
                    {
                        if (calculate_r_peptides_classification_data_template != null)
                        {
                            var template = calculate_r_peptides_classification_data_template.Select(a => new subsequence_classification_data.feature_info(a) {source = "", feature_value = 0}).ToList();

                            return template;
                        }
                        else
                        {
                            throw new Exception();
                        }
                    }


                    var aaComp_result = aaComp(engine, seq); // OAAC
                    //var aaDescriptors_result = aaDescriptors(seq); // variable length feature

                    //var lengthpep_result = lengthpep(seq);

                    var aIndex_result = aIndex(engine, seq);

                    var lag_start = 1;
                    var lag_end = 1;

                    var ph_start = 0;
                    var ph_end = 14;

                    var autoCorrelation_result_list = autoCorrelation_properties.SelectMany(b => Enumerable.Range(lag_start, (lag_end - lag_start) + 1).Select(c => (prop: b, lag: c, value: autoCorrelation(engine, seq, c, b, true))).ToList()).ToList();
                    var autoCovariance_result_list = autoCorrelation_properties.SelectMany(b => Enumerable.Range(lag_start, (lag_end - lag_start) + 1).Select(c => (prop: b, lag: c, value: autoCovariance(engine, seq, c, b, true))).ToList()).ToList();

                    var blosumIndices_result = blosumIndices(engine, seq);
                    blosumIndices_result.Add((-1, "Average", blosumIndices_result.Average(a => a.value)));

                    var boman_result = boman(engine, seq);
                    var charge_result_list = pKscales.SelectMany(a => Enumerable.Range(ph_start, (ph_end - ph_start) + 1).Select(b => (ph: b, scale: a, value: charge(engine, seq, b, a))).ToList()).ToList();
                    charge_result_list.AddRange(charge_result_list.GroupBy(a => a.scale).Select(a => (ph: -1, scale: $"{a.Key}_average", value: a.Select(b => b.value).Average())).ToList());

                    var crossCovariance_result_list = autoCorrelation_properties.SelectMany(a => autoCorrelation_properties.SelectMany(b => Enumerable.Range(lag_start, (lag_end - lag_start) + 1).Select(c => (prop1: a, prop2: b, lag: c, value: crossCovariance(engine, seq, c, a, b, true))).ToList()).ToList()).ToList();

                    var crucianiProperties_result = crucianiProperties(engine, seq); //3
                    var fasgaiVectors_result = fasgaiVectors(engine, seq); //6
                    var hmoment_result = hmoment(engine, seq);

                    var hydrophobicity_result_list = hydrophobicity_scales.Select(a => (scale: a, value: hydrophobicity(engine, seq, a))).ToList();


                    var instaIndex_result = instaIndex(engine, seq);
                    var kideraFactors_result = kideraFactors(engine, seq); //10

                    List<(int row, string Pwp, double H, double uH, double MembPosValue, string MembPos)> membpos_result = membpos(engine, seq);
                    membpos_result.Add((row: -1, Pwp: "", H: membpos_result.Count == 0 ? 0 : membpos_result.Average(a => a.H), uH: membpos_result.Count == 0 ? 0 : membpos_result.Average(a => a.uH), MembPosValue: membpos_result.Count == 0 ? 0 : membpos_result.Average(a => a.MembPosValue), MembPos: ""));


                    var mswhimScores_result = mswhimScores(engine, seq); //3
                    var mw_result = mw(engine, seq);
                    var pI_result = pI(engine, seq);
                    var protFP_result = protFP(engine, seq); //8
                    var stScales_result = stScales(engine, seq); //8
                    var tScales_result = tScales(engine, seq); //5
                    var vhseScales_result = vhseScales(engine, seq); //8
                    var zScales_result = zScales(engine, seq); //5


                    var features = new List<subsequence_classification_data.feature_info>();

                    var source = "";
                    var alpha_name = "Overall";

                    features.Add(new subsequence_classification_data.feature_info()
                    {
                        alphabet = alpha_name,
                        dimension = 1,
                        category = $"{nameof(r_peptides)}",
                        source = source,
                        @group = $"{nameof(r_peptides)}_{nameof(aIndex)}",
                        member = $"{nameof(r_peptides)}_{nameof(aIndex)}",
                        perspective = "default",
                        feature_value = descriptive_stats.fix_double(aIndex_result)
                    });

                    autoCorrelation_result_list.ForEach(autoCorrelation_result => features.Add(new subsequence_classification_data.feature_info()
                    {
                        alphabet = alpha_name,
                        dimension = 1,
                        category = $"{nameof(r_peptides)}",
                        source = source,
                        @group = $"{nameof(r_peptides)}_{nameof(autoCorrelation)}",
                        member = $"{nameof(r_peptides)}_{nameof(autoCorrelation)}_{autoCorrelation_result.prop}_lag{autoCorrelation_result.lag}",
                        perspective = "default",
                        feature_value = descriptive_stats.fix_double(autoCorrelation_result.value)
                    }));

                    autoCovariance_result_list.ForEach(autoCovariance_result => features.Add(new subsequence_classification_data.feature_info()
                    {
                        alphabet = alpha_name,
                        dimension = 1,
                        category = $"{nameof(r_peptides)}",
                        source = source,
                        @group = $"{nameof(r_peptides)}_{nameof(autoCovariance)}",
                        member = $"{nameof(r_peptides)}_{nameof(autoCovariance)}_{autoCovariance_result.prop}_lag{autoCovariance_result.lag}",
                        perspective = "default",
                        feature_value = descriptive_stats.fix_double(autoCovariance_result.value)
                    }));

                    features.AddRange(blosumIndices_result.Select(a => new subsequence_classification_data.feature_info()
                    {
                        alphabet = alpha_name,
                        dimension = 1,
                        category = $"{nameof(r_peptides)}",
                        source = source,
                        @group = $"{nameof(r_peptides)}_{nameof(blosumIndices)}",
                        member = $"{nameof(r_peptides)}_{nameof(blosumIndices)}_{a.index}_{a.name}",
                        perspective = "default",
                        feature_value = descriptive_stats.fix_double(a.value)
                    }).ToList());

                    features.Add(new subsequence_classification_data.feature_info()
                    {
                        alphabet = alpha_name,
                        dimension = 1,
                        category = $"{nameof(r_peptides)}",
                        source = source,
                        @group = $"{nameof(r_peptides)}_{nameof(boman)}",
                        member = $"{nameof(r_peptides)}_{nameof(boman)}",
                        perspective = "default",
                        feature_value = descriptive_stats.fix_double(boman_result)
                    });

                    charge_result_list.ForEach(charge_result => features.Add(new subsequence_classification_data.feature_info()
                    {
                        alphabet = alpha_name,
                        dimension = 1,
                        category = $"{nameof(r_peptides)}",
                        source = source,
                        @group = $"{nameof(r_peptides)}_{nameof(charge)}",
                        member = $"{nameof(r_peptides)}_{nameof(charge)}_{charge_result.scale}_ph{charge_result.ph}",
                        perspective = "default",
                        feature_value = descriptive_stats.fix_double(charge_result.value)
                    }));

                    crossCovariance_result_list.ForEach(crossCovariance_result => features.Add(new subsequence_classification_data.feature_info()
                    {
                        alphabet = alpha_name,
                        dimension = 1,
                        category = $"{nameof(r_peptides)}",
                        source = source,
                        @group = $"{nameof(r_peptides)}_{nameof(crossCovariance)}",
                        member = $"{nameof(r_peptides)}_{nameof(crossCovariance)}_{crossCovariance_result.prop1}_{crossCovariance_result.prop2}_lag{crossCovariance_result.lag}",
                        perspective = "default",
                        feature_value = descriptive_stats.fix_double(crossCovariance_result.value)
                    }));

                    features.AddRange(crucianiProperties_result.Select(a => new subsequence_classification_data.feature_info()
                    {
                        alphabet = alpha_name,
                        dimension = 1,
                        category = $"{nameof(r_peptides)}",
                        source = source,
                        @group = $"{nameof(r_peptides)}_{nameof(crucianiProperties)}",
                        member = $"{nameof(r_peptides)}_{nameof(crucianiProperties)}_{a.index}_{a.name}",
                        perspective = "default",
                        feature_value = descriptive_stats.fix_double(a.value)
                    }).ToList());

                    features.AddRange(fasgaiVectors_result.Select(a => new subsequence_classification_data.feature_info()
                    {
                        alphabet = alpha_name,
                        dimension = 1,
                        category = $"{nameof(r_peptides)}",
                        source = source,
                        @group = $"{nameof(r_peptides)}_{nameof(fasgaiVectors)}",
                        member = $"{nameof(r_peptides)}_{nameof(fasgaiVectors)}_{a.index}_{a.name}",
                        perspective = "default",
                        feature_value = descriptive_stats.fix_double(a.value)
                    }).ToList());

                    features.Add(new subsequence_classification_data.feature_info()
                    {
                        alphabet = alpha_name,
                        dimension = 1,
                        category = $"{nameof(r_peptides)}",
                        source = source,
                        @group = $"{nameof(r_peptides)}_{nameof(hmoment)}",
                        member = $"{nameof(r_peptides)}_{nameof(hmoment)}",
                        perspective = "default",
                        feature_value = descriptive_stats.fix_double(hmoment_result)
                    });

                    hydrophobicity_result_list.ForEach(hydrophobicity_result => features.Add(new subsequence_classification_data.feature_info()
                    {
                        alphabet = alpha_name,
                        dimension = 1,
                        category = $"{nameof(r_peptides)}",
                        source = source,
                        @group = $"{nameof(r_peptides)}_{nameof(hydrophobicity)}",
                        member = $"{nameof(r_peptides)}_{nameof(hydrophobicity)}_{hydrophobicity_result.scale}",
                        perspective = "default",
                        feature_value = descriptive_stats.fix_double(hydrophobicity_result.value)
                    }));

                    features.Add(new subsequence_classification_data.feature_info()
                    {
                        alphabet = alpha_name,
                        dimension = 1,
                        category = $"{nameof(r_peptides)}",
                        source = source,
                        @group = $"{nameof(r_peptides)}_{nameof(instaIndex)}",
                        member = $"{nameof(r_peptides)}_{nameof(instaIndex)}",
                        perspective = "default",
                        feature_value = descriptive_stats.fix_double(instaIndex_result)
                    });

                    features.AddRange(kideraFactors_result.Select(a => new subsequence_classification_data.feature_info()
                    {
                        alphabet = alpha_name,
                        dimension = 1,
                        category = $"{nameof(r_peptides)}",
                        source = source,
                        @group = $"{nameof(r_peptides)}_{nameof(kideraFactors)}",
                        member = $"{nameof(r_peptides)}_{nameof(kideraFactors)}_{a.index}_{a.name}",
                        perspective = "default",
                        feature_value = descriptive_stats.fix_double(a.value)
                    }).ToList());


                    //features.AddRange(membpos_result_av.Select(a => new subsequence_classification_data.feature_info()
                    //{
                    //    alphabet = alpha_name,
                    //    dimension = 1,
                    //    category = $"{nameof(r_peptides)}",
                    //    source = source,
                    //    @group = $"{nameof(r_peptides)}_{nameof(membpos)}",
                    //    member = $"{nameof(r_peptides)}_{nameof(membpos)}_{a.name}",
                    //    perspective = "default",
                    //    feature_value = descriptive_stats.fix_double(a.value)
                    //}).ToList());

                    features.AddRange(mswhimScores_result.Select(a => new subsequence_classification_data.feature_info()
                    {
                        alphabet = alpha_name,
                        dimension = 1,
                        category = $"{nameof(r_peptides)}",
                        source = source,
                        @group = $"{nameof(r_peptides)}_{nameof(mswhimScores)}",
                        member = $"{nameof(r_peptides)}_{nameof(mswhimScores)}_{a.index}_{a.name}",
                        perspective = "default",
                        feature_value = descriptive_stats.fix_double(a.value)
                    }).ToList());

                    features.Add(new subsequence_classification_data.feature_info()
                    {
                        alphabet = alpha_name,
                        dimension = 1,
                        category = $"{nameof(r_peptides)}",
                        source = source,
                        @group = $"{nameof(r_peptides)}_{nameof(mw)}",
                        member = $"{nameof(r_peptides)}_{nameof(mw)}",
                        perspective = "default",
                        feature_value = descriptive_stats.fix_double(mw_result)
                    });

                    features.Add(new subsequence_classification_data.feature_info()
                    {
                        alphabet = alpha_name,
                        dimension = 1,
                        category = $"{nameof(r_peptides)}",
                        source = source,
                        @group = $"{nameof(r_peptides)}_{nameof(pI)}",
                        member = $"{nameof(r_peptides)}_{nameof(pI)}",
                        perspective = "default",
                        feature_value = descriptive_stats.fix_double(pI_result)
                    });

                    features.AddRange(protFP_result.Select(a => new subsequence_classification_data.feature_info()
                    {
                        alphabet = alpha_name,
                        dimension = 1,
                        category = $"{nameof(r_peptides)}",
                        source = source,
                        @group = $"{nameof(r_peptides)}_{nameof(protFP)}",
                        member = $"{nameof(r_peptides)}_{nameof(protFP)}_{a.index}_{a.name}",
                        perspective = "default",
                        feature_value = descriptive_stats.fix_double(a.value)
                    }).ToList());

                    features.AddRange(stScales_result.Select(a => new subsequence_classification_data.feature_info()
                    {
                        alphabet = alpha_name,
                        dimension = 1,
                        category = $"{nameof(r_peptides)}",
                        source = source,
                        @group = $"{nameof(r_peptides)}_{nameof(stScales)}",
                        member = $"{nameof(r_peptides)}_{nameof(stScales)}_{a.index}_{a.name}",
                        perspective = "default",
                        feature_value = descriptive_stats.fix_double(a.value)
                    }).ToList());

                    features.AddRange(tScales_result.Select(a => new subsequence_classification_data.feature_info()
                    {
                        alphabet = alpha_name,
                        dimension = 1,
                        category = $"{nameof(r_peptides)}",
                        source = source,
                        @group = $"{nameof(r_peptides)}_{nameof(tScales)}",
                        member = $"{nameof(r_peptides)}_{nameof(tScales)}_{a.index}_{a.name}",
                        perspective = "default",
                        feature_value = descriptive_stats.fix_double(a.value)
                    }).ToList());

                    features.AddRange(vhseScales_result.Select(a => new subsequence_classification_data.feature_info()
                    {
                        alphabet = alpha_name,
                        dimension = 1,
                        category = $"{nameof(r_peptides)}",
                        source = source,
                        @group = $"{nameof(r_peptides)}_{nameof(vhseScales)}",
                        member = $"{nameof(r_peptides)}_{nameof(vhseScales)}_{a.index}_{a.name}",
                        perspective = "default",
                        feature_value = descriptive_stats.fix_double(a.value)
                    }).ToList());

                    features.AddRange(zScales_result.Select(a => new subsequence_classification_data.feature_info()
                    {
                        alphabet = alpha_name,
                        dimension = 1,
                        category = $"{nameof(r_peptides)}",
                        source = source,
                        @group = $"{nameof(r_peptides)}_{nameof(zScales)}",
                        member = $"{nameof(r_peptides)}_{nameof(zScales)}_{a.index}_{a.name}",
                        perspective = "default",
                        feature_value = descriptive_stats.fix_double(a.value)
                    }).ToList());




                    //Console.WriteLine($"{nameof(aaComp_result)} = {string.Join(",", aaComp_result)}");
                    //Console.WriteLine($"{nameof(aaDescriptors_result)} = {string.Join(",", aaDescriptors_result)}");

                    //Console.WriteLine($"{nameof(aIndex_result)} = {string.Join(",", aIndex_result)}");
                    //Console.WriteLine($"{nameof(autoCorrelation_result)} = {string.Join(",", autoCorrelation_result)}");
                    //Console.WriteLine($"{nameof(autoCovariance_result)} = {string.Join(",", autoCovariance_result)}");
                    //Console.WriteLine($"{nameof(blosumIndices_result)} = {string.Join(",", blosumIndices_result)}");
                    //Console.WriteLine($"{nameof(boman_result)} = {string.Join(",", boman_result)}");
                    //Console.WriteLine($"{nameof(charge_result)} = {string.Join(",", charge_result)}");
                    //Console.WriteLine($"{nameof(crossCovariance_result)} = {string.Join(",", crossCovariance_result)}");
                    //Console.WriteLine($"{nameof(crucianiProperties_result)} = {string.Join(",", crucianiProperties_result)}");
                    //Console.WriteLine($"{nameof(fasgaiVectors_result)} = {string.Join(",", fasgaiVectors_result)}");
                    //Console.WriteLine($"{nameof(hmoment_result)} = {string.Join(",", hmoment_result)}");
                    //Console.WriteLine($"{nameof(hydrophobicity_result)} = {string.Join(",", hydrophobicity_result)}");
                    //Console.WriteLine($"{nameof(instaIndex_result)} = {string.Join(",", instaIndex_result)}");
                    //Console.WriteLine($"{nameof(kideraFactors_result)} = {string.Join(",", kideraFactors_result)}");
                    //Console.WriteLine($"{nameof(lengthpep_result)} = {string.Join(",", lengthpep_result)}");
                    //Console.WriteLine($"{nameof(membpos_result)} = {string.Join(",", membpos_result)}");
                    //Console.WriteLine($"{nameof(mswhimScores_result)} = {string.Join(",", mswhimScores_result)}");
                    //Console.WriteLine($"{nameof(mw_result)} = {string.Join(",", mw_result)}");
                    //Console.WriteLine($"{nameof(pI_result)} = {string.Join(",", pI_result)}");
                    //Console.WriteLine($"{nameof(protFP_result)} = {string.Join(",", protFP_result)}");
                    //Console.WriteLine($"{nameof(stScales_result)} = {string.Join(",", stScales_result)}");
                    //Console.WriteLine($"{nameof(tScales_result)} = {string.Join(",", tScales_result)}");
                    //Console.WriteLine($"{nameof(vhseScales_result)} = {string.Join(",", vhseScales_result)}");
                    //Console.WriteLine($"{nameof(zScales_result)} = {string.Join(",", zScales_result)}");
                    //Console.WriteLine("");

                    if (calculate_r_peptides_classification_data_template == null)
                    {
                        var template = features.Select(a => new subsequence_classification_data.feature_info(a) {source = "", feature_value = 0}).ToList();
                        calculate_r_peptides_classification_data_template = template;
                    }



                    return features;
                }
            }
        }

        /// <summary>
        /// This function calculates the amount of amino acids of a particular class and classified as: Tiny, Small, Aliphatic, Aromatic, Non-polar, Polar, Charged, Basic and Acidic based on their size and R-groups using same function implemented in EMBOSS 'pepstat'. The output is a matrix with the number and percentage of amino acids of a particular class
        ///
        /// This function was originally written by Alan Bleasby (ajb@ebi.ac.uk) for the EMBOSS package. Further information: http://emboss.sourceforge.net/apps/cvs/emboss/apps/pepstats.html
        ///
        /// Rice, Peter, Ian Longden, and Alan Bleasby. "EMBOSS: the European molecular biology open software suite." Trends in genetics 16.6 (2000): 276-277.
        /// </summary>
        /// <param name="seq">An amino-acid sequence</param>
        /// <returns>
        /// The output is a matrix with the number and percentage of amino acids of a particular class
        /// 
        /// Tiny (A + C + G + S + T)
        /// Small (A + B + C + D + G + N + P + S + T + V)
        /// Aliphatic (A + I + L + V)
        /// Aromatic (F + H + W + Y)
        /// Non-polar (A + C + F + G + I + L + M + P + V + W + Y)
        /// Polar (D + E + H + K + N + Q + R + S + T + Z)
        /// Charged (B + D + E + H + K + R + Z)
        /// Basic (H + K + R)
        /// Acidic (B + D + E + Z)
        /// </returns>
        public static List<(int row, string rowname, double count, double value)> aaComp(REngine engine, string seq)
        {
            // https://rdrr.io/cran/Peptides/man/aaComp.html

            lock (engine_lock)
            {
                var f = nameof(aaComp);
                var k = Key;
                var v = $"{f}_v{k}";

                var evaluate = engine.Evaluate($"{v} <- {f}(seq = \"{seq}\")");
                var rownames = engine.Evaluate($"rownames({v}[[1]])").AsCharacter();
                var colnames = engine.Evaluate($"colnames({v}[[1]])").AsCharacter();
                var values = engine.Evaluate($"{v}[[1]]").AsNumericMatrix();
                var rm = engine.Evaluate($"rm({v})");


                var list = new List<(int row, string rowname, double count, double value)>();


                for (var r = 0; r < rownames.Length; r++)

                {
                    var x = (r, rownames[r], values[r, 0], values[r, 1]);
                    list.Add(x);
                }



                //Tiny(A + C + G + S + T)
                //Small(A + B + C + D + G + N + P + S + T + V)
                //Aliphatic(A + I + L + V)
                //Aromatic(F + H + W + Y)
                //Non - polar(A + C + F + G + I + L + M + P + V + W + Y)
                //Polar(D + E + H + K + N + Q + R + S + T + Z)
                //Charged(B + D + E + H + K + R + Z)
                //Basic(H + K + R)
                //Acidic(B + D + E + Z)

                return list;
            }
        }

        /// <summary>
        /// 
        /// aaDescriptors: Compute 66 descriptors for each amino acid of a protein...
        /// In Peptides: Calculate Indices and Theoretical Physicochemical Properties of Protein Sequences
        /// Description Usage Arguments Value Examples
        /// 
        /// View source: R/aaDescriptors.R
        /// 
        /// Description
        /// The function return 66 amino acid descriptors for the 20 natural amino acids. Available descriptors are:
        /// 
        /// crucianiProperties: Cruciani, G., Baroni, M., Carosati, E., Clementi, M., Valigi, R., and Clementi, S. (2004) Peptide studies by means of principal properties of amino acids derived from MIF descriptors. J. Chemom. 18, 146-155.,
        /// kideraFactors: Kidera, A., Konishi, Y., Oka, M., Ooi, T., &amp; Scheraga, H. A. (1985). Statistical analysis of the physical properties of the 20 naturally occurring amino acids. Journal of Protein Chemistry, 4(1), 23-55.,
        /// zScales: Sandberg M, Eriksson L, Jonsson J, Sjostrom M, Wold S: New chemical descriptors relevant for the design of biologically active peptides. A multivariate characterization of 87 amino acids. J Med Chem 1998, 41:2481-2491.,
        /// FASGAI: Liang, G., &amp; Li, Z. (2007). Factor analysis scale of generalized amino acid information as the source of a new set of descriptors for elucidating the structure and activity relationships of cationic antimicrobial peptides. Molecular Informatics, 26(6), 754-763.,
        /// tScales: Tian F, Zhou P, Li Z: T-scale as a novel vector of topological descriptors for amino acids and its application in QSARs of peptides. J Mol Struct. 2007, 830: 106-115. 10.1016/j.molstruc.2006.07.004.,
        /// VHSE: VHSE-scales (principal components score Vectors of Hydrophobic, Steric, and Electronic properties), is derived from principal components analysis (PCA) on independent families of 18 hydrophobic properties, 17 steric properties, and 15 electronic properties, respectively, which are included in total 50 physicochemical variables of 20 coded amino acids.,
        /// protFP: van Westen, G. J., Swier, R. F., Wegner, J. K., IJzerman, A. P., van Vlijmen, H. W., &amp; Bender, A. (2013). Benchmarking of protein descriptor sets in proteochemometric modeling (part 1): comparative study of 13 amino acid descriptor sets. Journal of cheminformatics, 5(1), 41.,
        /// stScales: Yang, L., Shu, M., Ma, K., Mei, H., Jiang, Y., &amp; Li, Z. (2010). ST-scale as a novel amino acid descriptor and its application in QSAM of peptides and analogues. Amino acids, 38(3), 805-816.,
        /// BLOSUM: Georgiev, A. G. (2009). Interpretable numerical descriptors of amino acid space. Journal of Computational Biology, 16(5), 703-723.,
        /// MSWHIM: Zaliani, A., &amp; Gancia, E. (1999). MS-WHIM scores for amino acids: a new 3D-description for peptide QSAR and QSPR studies. Journal of chemical information and computer sciences, 39(3), 525-533.
        ///
        /// https://rdrr.io/cran/Peptides/man/aaDescriptors.html
        /// </summary>
        /// <param name="seq">An amino-acids sequence. If multiple sequences are given all of them must have the same length (gap symbols are allowed.)</param>
        /// <returns>a matrix with 66 amino acid descriptors for each aminoacid in a protein sequence.</returns>
        public static List<(int row, string rowname, int col, string colname, double value)> aaDescriptors(REngine engine, string seq)
        {
            lock (engine_lock)
            {
                var f = nameof(aaDescriptors);
                var k = Key;
                var v = $"{f}_v{k}";

                var evaluate = engine.Evaluate($"{v} <- {f}(seq = \"{seq}\")");

                var colnames = engine.Evaluate($"colnames({v})").AsCharacter();
                var values = engine.Evaluate($"{v}").AsNumeric();
                var rm = engine.Evaluate($"rm({v})");

                var list = new List<(int row, string rowname, int col, string colname, double value)>();

                var total_rows = seq.Length;
                var total_cols = values.Length / seq.Length;

                var i = 0;

                for (var r = 0; r < total_rows; r++)
                {
                    for (var c = 0; c < total_cols; c++)

                    {



                        var x = (r, "", c, colnames[i], values[i]);

                        list.Add(x);

                        i++;
                    }
                }

                return list;
            }
        }

        /// <summary>
        /// This function calculates the Ikai(1980) aliphatic index of a protein.The aindex is defined as the relative volume occupied by aliphatic side chains(Alanine, Valine, Isoleucine, and Leucine). It may be regarded as a positive factor for the increase of thermostability of globular proteins.
        ///
        /// Ikai (1980). Thermostability and aliphatic index of globular proteins. Journal of Biochemistry, 88(6), 1895-1898.
        /// </summary>
        /// <param name="seq">An amino-acids sequence</param>
        /// <returns>The computed aliphatic index for a given amino-acids sequence</returns>
        public static double aIndex(REngine engine, string seq)
        {
            lock (engine_lock)
            {
                var f = nameof(aIndex);
                var k = Key;
                var v = $"{f}_v{k}";

                var evaluate = engine.Evaluate($"{v} <- {f}(seq = \"{seq}\")");

                var values = engine.Evaluate($"{v}[[1]]").AsNumeric();
                var rm = engine.Evaluate($"rm({v})");

                var result = values.First();



                return result;
            }
        }

        public static string[] autoCorrelation_properties = new string[]
            {
//"$BLOSUM",
"$BLOSUM$BLOSUM1",
"$BLOSUM$BLOSUM10",
"$BLOSUM$BLOSUM2",
"$BLOSUM$BLOSUM3",
"$BLOSUM$BLOSUM4",
"$BLOSUM$BLOSUM5",
"$BLOSUM$BLOSUM6",
"$BLOSUM$BLOSUM7",
"$BLOSUM$BLOSUM8",
"$BLOSUM$BLOSUM9",
//"$crucianiProperties",
"$crucianiProperties$PP1",
"$crucianiProperties$PP2",
"$crucianiProperties$PP3",
//"$FASGAI",
"$FASGAI$F1",
"$FASGAI$F2",
"$FASGAI$F3",
"$FASGAI$F4",
"$FASGAI$F5",
"$FASGAI$F6",
//"$Hydrophobicity",
"$Hydrophobicity$Aboderin",
"$Hydrophobicity$AbrahamLeo",
"$Hydrophobicity$Argos",
"$Hydrophobicity$BlackMould",
"$Hydrophobicity$BullBreese",
"$Hydrophobicity$Casari",
"$Hydrophobicity$Chothia",
"$Hydrophobicity$Cid",
"$Hydrophobicity$Cowan3.4",
"$Hydrophobicity$Cowan7.5",
"$Hydrophobicity$Eisenberg",
"$Hydrophobicity$Engelman",
"$Hydrophobicity$Fasman",
"$Hydrophobicity$Fauchere",
"$Hydrophobicity$Goldsack",
"$Hydrophobicity$Guy",
"$Hydrophobicity$HoppWoods",
"$Hydrophobicity$interfaceScale_pH2",
"$Hydrophobicity$interfaceScale_pH8",
"$Hydrophobicity$Janin",
"$Hydrophobicity$Jones",
"$Hydrophobicity$Juretic",
"$Hydrophobicity$Kidera",
"$Hydrophobicity$Kuhn",
"$Hydrophobicity$KyteDoolittle",
"$Hydrophobicity$Levitt",
"$Hydrophobicity$Manavalan",
"$Hydrophobicity$Miyazawa",
"$Hydrophobicity$octanolScale_pH2",
"$Hydrophobicity$octanolScale_pH8",
"$Hydrophobicity$oiScale_pH2",
"$Hydrophobicity$oiScale_pH8",
"$Hydrophobicity$Parker",
"$Hydrophobicity$Ponnuswamy",
"$Hydrophobicity$Prabhakaran",
"$Hydrophobicity$Rao",
"$Hydrophobicity$Rose",
"$Hydrophobicity$Roseman",
"$Hydrophobicity$Sweet",
"$Hydrophobicity$Tanford",
"$Hydrophobicity$Welling",
"$Hydrophobicity$Wilson",
"$Hydrophobicity$Wolfenden",
"$Hydrophobicity$Zimmerman",
//"$kideraFactors",
"$kideraFactors$KF1",
"$kideraFactors$KF10",
"$kideraFactors$KF2",
"$kideraFactors$KF3",
"$kideraFactors$KF4",
"$kideraFactors$KF5",
"$kideraFactors$KF6",
"$kideraFactors$KF7",
"$kideraFactors$KF8",
"$kideraFactors$KF9",
//"$MSWHIM",
"$MSWHIM$MSWHIM1",
"$MSWHIM$MSWHIM2",
"$MSWHIM$MSWHIM3",
//"$pK",
"$pK$Bjellqvist",
"$pK$Dawson",
"$pK$EMBOSS",
"$pK$Lehninger",
"$pK$Murray",
"$pK$Rodwell",
"$pK$Sillero",
"$pK$Solomon",
"$pK$Stryer",
//"$ProtFP",
"$ProtFP$ProtFP1",
"$ProtFP$ProtFP2",
"$ProtFP$ProtFP3",
"$ProtFP$ProtFP4",
"$ProtFP$ProtFP5",
"$ProtFP$ProtFP6",
"$ProtFP$ProtFP7",
"$ProtFP$ProtFP8",
//"$stScales",
"$stScales$ST1",
"$stScales$ST2",
"$stScales$ST3",
"$stScales$ST4",
"$stScales$ST5",
"$stScales$ST6",
"$stScales$ST7",
"$stScales$ST8",
//"$tScales",
"$tScales$T1",
"$tScales$T2",
"$tScales$T3",
"$tScales$T4",
"$tScales$T5",
//"$VHSE",
"$VHSE$VHSE1",
"$VHSE$VHSE2",
"$VHSE$VHSE3",
"$VHSE$VHSE4",
"$VHSE$VHSE5",
"$VHSE$VHSE6",
"$VHSE$VHSE7",
"$VHSE$VHSE8",
//"$zScales",
"$zScales$Z1",
"$zScales$Z2",
"$zScales$Z3",
"$zScales$Z4",
"$zScales$Z5"
};

        /// <summary>
        /// This function computes the Cruciani et al (2004) auto-correlation index. The autoCorrelation index is calculated for a lag 'd' using a descriptor 'f' (centred) over a sequence of length 'L'.
        ///
        /// Cruciani, G., Baroni, M., Carosati, E., Clementi, M., Valigi, R., and Clementi, S. (2004) Peptide studies by means of principal properties of amino acids derived from MIF descriptors. J. Chemom. 18, 146-155.
        ///
        /// https://rdrr.io/cran/Peptides/man/autoCorrelation.html
        /// </summary>
        /// <param name="seq">An amino-acids sequence</param>
        /// <param name="lag">A value for a lag, the max value is equal to the length of shortest peptide minus one.</param>
        /// <param name="property">A property to use as value to be correlated.</param>
        /// <param name="center">A logical value TRUE or FALSE if the property must be centered.</param>
        /// <returns>The computed auto-correlation index for a given amino-acids sequence</returns>
        public static double autoCorrelation(REngine engine, string seq, int lag, string property, bool center = true)
        {
            lock (engine_lock)
            {
                if (!autoCorrelation_properties.Contains(property)) throw new Exception();

                if (lag >= seq.Length - 1) return 0;

                var f = nameof(autoCorrelation);
                var k = Key;
                var v = $"{f}_v{k}";

                var aadata = engine.Evaluate("if (!exists(\"AAdata\")) data(AAdata)");
                var evaluate = engine.Evaluate($"{v} <- {f}(sequence = \"{seq}\", lag = {lag}, property = AAdata{property}, center = {center.ToString().ToUpperInvariant()})");

                var values = engine.Evaluate($"{v}[[1]]").AsNumeric();
                var rm = engine.Evaluate($"rm({v})");

                var result = values.First();


                return result;
            }
        }

        /// <summary>
        /// This function computes the Cruciani et al (2004) auto-corvariance index. The autoCovariance index is calculated for a lag 'd' using a descriptor 'f' (centred) over a sequence of length 'L'.
        ///
        /// Cruciani, G., Baroni, M., Carosati, E., Clementi, M., Valigi, R., and Clementi, S. (2004) Peptide studies by means of principal properties of amino acids derived from MIF descriptors. J. Chemom. 18, 146-155.
        /// </summary>
        /// <param name="seq">An amino-acids sequence</param>
        /// <param name="lag">A value for a lag, the max value is equal to the length of the shortest peptide minus one.</param>
        /// <param name="property">A property to use as value to evaluate the covariance.</param>
        /// <param name="center">A logical value TRUE or FALSE if the property must be centered.</param>
        /// <returns>The computed auto-covariance index for a given amino-acids sequence</returns>
        public static double autoCovariance(REngine engine, string seq, int lag, string property, bool center = true)
        {
            lock (engine_lock)
            {
                if (!autoCorrelation_properties.Contains(property)) throw new Exception();
                if (lag >= seq.Length - 1) return 0;


                var f = nameof(autoCovariance);
                var k = Key;
                var v = $"{f}_v{k}";

                var aadata = engine.Evaluate("if (!exists(\"AAdata\")) data(AAdata)");
                var evaluate = engine.Evaluate($"{v} <- {f}(sequence = \"{seq}\", lag = {lag}, property = AAdata{property}, center = {center.ToString().ToUpperInvariant()})");

                var values = engine.Evaluate($"{v}[[1]]").AsNumeric();
                var rm = engine.Evaluate($"rm({v})");

                var result = values.First();


                return result;
            }
        }

        /// <summary>
        /// BLOSUM indices were derived of physicochemical properties that have been subjected to a VARIMAX analyses and an alignment matrix of the 20 natural AAs using the BLOSUM62 matrix.
        ///
        /// Georgiev, A. G. (2009). Interpretable numerical descriptors of amino acid space. Journal of Computational Biology, 16(5), 703-723.
        /// </summary>
        /// <param name="seq">An amino-acids sequence</param>
        /// <returns>The computed average of BLOSUM indices of all the amino acids in the corresponding peptide sequence.</returns>
        public static List<(int index, string name, double value)> blosumIndices(REngine engine, string seq)
        {
            lock (engine_lock)
            {
                var f = nameof(blosumIndices);
                var k = Key;
                var v = $"{f}_v{k}";
                var ai = "[[1]]";

                var evaluate = engine.Evaluate($"{v} <- {f}(seq = \"{seq}\")");
                var names = engine.Evaluate($"names({v}{ai})").AsCharacter();
                var dimnames = engine.Evaluate($"dimnames({v}{ai})");
                var rownames = engine.Evaluate($"rownames({v}{ai})");
                var colnames = engine.Evaluate($"colnames({v}{ai})");
                var values = engine.Evaluate($"{v}{ai}").AsNumeric();
                var rm = engine.Evaluate($"rm({v})");

                var result = new List<(int index, string name, double value)>();

                for (var i = 0; i < names.Length; i++)
                {
                    var x = (i, names[i], values[i]);
                    result.Add(x);
                }

                return result;
            }
        }

        /// <summary>
        /// This function computes the potential protein interaction index proposed by Boman (2003) based in the amino acid sequence of a protein. The index is equal to the sum of the solubility values for all residues in a sequence, it might give an overall estimate of the potential of a peptide to bind to membranes or other proteins as receptors, to normalize it is divided by the number of residues. A protein have high binding potential if the index value is higher than 2.48.
        ///
        /// The potential protein interaction index was proposed by Boman (2003) as an easy way to differentiate the action mechanism of hormones (protein-protein) and antimicrobial peptides (protein-membrane) through this index. This function predicts the potential peptide interaction with another protein.
        /// 
        /// Boman, H. G. (2003). Antibacterial peptides: basic facts and emerging concepts. Journal of Internal Medicine, 254(3), 197-215.
        /// </summary>
        /// <param name="seq">An amino-acid sequence</param>
        /// <returns>The computed potential protein-protein interaction for a given amino-acids sequence</returns>
        public static double boman(REngine engine, string seq)
        {
            lock (engine_lock)
            {
                var f = nameof(boman);
                var k = Key;
                var v = $"{f}_v{k}";

                var evaluate = engine.Evaluate($"{v} <- {f}(seq = \"{seq}\")");

                var values = engine.Evaluate($"{v}[[1]]").AsNumeric();
                var rm = engine.Evaluate($"rm({v})");

                var result = values.First();


                return result;
            }
        }

        public static string[] pKscales = new string[]
        {
            "Bjellqvist", "Dawson", "EMBOSS", "Lehninger", "Murray", "Rodwell", "Sillero", "Solomon", "Stryer"
        };


        /// <summary>
        /// This function computes the net charge of a protein sequence based on the Henderson-Hasselbalch equation described by Moore, D. S. (1985). The net charge can be calculated at defined pH using one of the 9 pKa scales availables: Bjellqvist, Dawson, EMBOSS, Lehninger, Murray, Rodwell, Sillero, Solomon or Stryer.
        ///
        /// Kiraga, J. (2008) Analysis and computer simulations of variability of isoelectric point of proteins in the proteomes. PhD thesis, University of Wroclaw, Poland.
        /// 
        /// Bjellqvist, B., Hughes, G.J., Pasquali, Ch., Paquet, N., Ravier, F., Sanchez, J.Ch., Frutige,r S., Hochstrasser D. (1993) The focusing positions of polypeptides in immobilized pH gradients can be predicted from their amino acid sequences. Electrophoresis, 14:1023-1031.
        /// 
        /// Dawson, R. M. C.; Elliot, D. C.; Elliot, W. H.; Jones, K. M. Data for biochemical research. Oxford University Press, 1989; p. 592.
        /// 
        /// EMBOSS data are from http://emboss.sourceforge.net/apps/release/5.0/emboss/apps/iep.html.
        /// 
        /// Nelson, D. L.; Cox, M. M. Lehninger Principles of Biochemistry, Fourth Edition; W. H. Freeman, 2004; p. 1100.
        /// 
        /// Murray, R.K., Granner, D.K., Rodwell, V.W. (2006) Harper's illustrated Biochemistry. 27th edition. Published by The McGraw-Hill Companies.
        /// 
        /// Rodwell, J. Heterogeneity of component bands in isoelectric focusing patterns. Analytical Biochemistry, 1982, 119 (2), 440-449.
        /// 
        /// Sillero, A., Maldonado, A. (2006) Isoelectric point determination of proteins and other macromolecules: oscillating method. Comput Biol Med., 36:157-166.
        /// 
        /// Solomon, T.W.G. (1998) Fundamentals of Organic Chemistry, 5th edition. Published by Wiley.
        /// 
        /// Stryer L. (1999) Biochemia. czwarta edycja. Wydawnictwo Naukowe PWN.
        /// </summary>
        /// <param name="seq">An amino-acids sequence</param>
        /// <param name="pH">A pH value</param>
        /// <param name="pKscale">A character string specifying the pKa scale to be used; must be one of "Bjellqvist", "Dawson", "EMBOSS", "Lehninger", "Murray", "Rodwell", "Sillero", "Solomon" or "Stryer"</param>
        /// <returns></returns>
        public static double charge(REngine engine, string seq, double pH = 7, string pKscale = "Lehninger")
        {
            lock (engine_lock)
            {
                if (!pKscales.Contains(pKscale)) throw new Exception();

                var f = nameof(charge);
                var k = Key;
                var v = $"{f}_v{k}";

                var evaluate = engine.Evaluate($"{v} <- {f}(seq = \"{seq}\", pH = {pH}, pKscale = \"{pKscale}\")");
                var values = engine.Evaluate($"{v}[[1]]").AsNumeric();
                var rm = engine.Evaluate($"rm({v})");

                var result = values.First();


                return result;
            }
        }

        /// <summary>
        /// This function computes the Cruciani et al (2004) cross-covariance index. The lagged crossCovariance index is calculated for a lag 'd' using two descriptors 'f1' and 'f2' (centred) over a sequence of length 'L'
        ///
        /// Cruciani, G., Baroni, M., Carosati, E., Clementi, M., Valigi, R., and Clementi, S. (2004) Peptide studies by means of principal properties of amino acids derived from MIF descriptors. J. Chemom. 18, 146-155.
        /// </summary>
        /// <param name="seq">An amino-acids sequence</param>
        /// <param name="lag">A value for a lag, the max value is equal to the length of the shortest peptide minus one</param>
        /// <param name="property1">A property to use as value to evaluate the cross-covariance.</param>
        /// <param name="property2">A property to use as value to evaluate the cross-covariance.</param>
        /// <param name="center">A logical value TRUE or FALSE if the property must be centered.</param>
        /// <returns>The computed cross-covariance index for a given amino-acids sequence</returns>
        public static double crossCovariance(REngine engine, string seq, int lag, string property1, string property2, bool center = true)
        {
            lock (engine_lock)
            {
                if (!autoCorrelation_properties.Contains(property1)) throw new Exception();
                if (!autoCorrelation_properties.Contains(property2)) throw new Exception();

                if (lag >= seq.Length - 1) return 0;


                var f = nameof(crossCovariance);
                var k = Key;
                var v = $"{f}_v{k}";

                var aadata = engine.Evaluate("if (!exists(\"AAdata\")) data(AAdata)");
                var evaluate = engine.Evaluate($"{v} <- {f}(sequence = \"{seq}\", lag = {lag}, property1 = AAdata{property1}, property2 = AAdata{property2}, center = {center.ToString().ToUpperInvariant()})");

                var values = engine.Evaluate($"{v}[[1]]").AsNumeric();
                var rm = engine.Evaluate($"rm({v})");

                var result = values.First();


                return result;
            }
        }

        /// <summary>
        /// This function calculates the Cruciani properties of an amino-acids sequence using the scaled principal component scores that summarize a broad set of descriptors calculated based on the interaction of each amino acid residue with several chemical groups(or "probes"), such as charged ions, methyl, hydroxyl groups, and so forth.
        ///
        /// Cruciani, G., Baroni, M., Carosati, E., Clementi, M., Valigi, R., and Clementi, S. (2004) Peptide studies by means of principal properties of amino acids derived from MIF descriptors. J. Chemom. 18, 146-155.
        /// </summary>
        /// <param name="seq">An amino-acids sequence</param>
        /// <returns>
        /// The computed average of Cruciani properties of all the amino acids in the corresponding peptide sequence. Each PP represent an amino-acid property as follows:
        /// PP1: Polarity,
        /// PP2: Hydrophobicity,
        /// PP3: H-bonding
        /// </returns>
        public static List<(int index, string name, string description, double value)> crucianiProperties(REngine engine, string seq)
        {
            lock (engine_lock)
            {
                var f = nameof(crucianiProperties);
                var k = Key;
                var v = $"{f}_v{k}";
                var ai = "[[1]]";

                var evaluate = engine.Evaluate($"{v} <- {f}(seq = \"{seq}\")");
                var names = engine.Evaluate($"names({v}{ai})").AsCharacter();
                var dimnames = engine.Evaluate($"dimnames({v}{ai})");
                var rownames = engine.Evaluate($"rownames({v}{ai})");
                var colnames = engine.Evaluate($"colnames({v}{ai})");
                var values = engine.Evaluate($"{v}{ai}").AsNumeric();
                var rm = engine.Evaluate($"rm({v})");

                var result = new List<(int index, string name, string description, double value)>();

                for (var i = 0; i < names.Length; i++)
                {
                    var name = names[i];
                    var desc = "";

                    if (name == "PP1") desc = "Polarity";
                    else if (name == "PP2") desc = "Hydrophobicity";
                    else if (name == "PP3") desc = "H-bonding";

                    var x = (i, names[i], desc, values[i]);
                    result.Add(x);
                }

                return result;
            }
        }

        /// <summary>
        /// The FASGAI vectors (Factor Analysis Scales of Generalized Amino Acid Information) is a set of amino acid descriptors, that reflects hydrophobicity, alpha and turn propensities, bulky properties, compositional characteristics, local flexibility, and electronic properties, that can be utilized to represent the sequence structural features of peptides or protein motifs.
        ///
        /// Liang, G., & Li, Z. (2007). Factor analysis scale of generalized amino acid information as the source of a new set of descriptors for elucidating the structure and activity relationships of cationic antimicrobial peptides. Molecular Informatics, 26(6), 754-763.
        /// </summary>
        /// <param name="seq">An amino-acids sequence</param>
        /// <returns>
        /// The computed average of FASGAI factors of all the amino acids in the corresponding peptide sequence. Each factor represent an amino-acid property as follows:
        /// 
        /// F1: Hydrophobicity index,
        /// 
        /// F2: Alpha and turn propensities,
        /// 
        /// F3: Bulky properties,
        /// 
        /// F4: Compositional characteristic index,
        /// 
        /// F5: Local flexibility,
        /// 
        /// F6: Electronic properties
        /// </returns>

        public static List<(int index, string name, string description, double value)> fasgaiVectors(REngine engine, string seq)
        {
            lock (engine_lock)
            {
                var f = nameof(fasgaiVectors);
                var k = Key;
                var v = $"{f}_v{k}";
                var ai = "[[1]]";

                var evaluate = engine.Evaluate($"{v} <- {f}(seq = \"{seq}\")");
                var names = engine.Evaluate($"names({v}{ai})").AsCharacter();
                var dimnames = engine.Evaluate($"dimnames({v}{ai})");
                var rownames = engine.Evaluate($"rownames({v}{ai})");
                var colnames = engine.Evaluate($"colnames({v}{ai})");
                var values = engine.Evaluate($"{v}{ai}").AsNumeric();
                var rm = engine.Evaluate($"rm({v})");

                var result = new List<(int index, string name, string description, double value)>();

                for (var i = 0; i < names.Length; i++)
                {
                    var name = names[i];
                    var desc = "";

                    if (name == "F1") desc = "Hydrophobicity index";
                    else if (name == "F2") desc = "Alpha and turn propensities";
                    else if (name == "F3") desc = "Bulky properties";
                    else if (name == "F4") desc = "Compositional characteristic index";
                    else if (name == "F5") desc = "Local flexibility";
                    else if (name == "F6") desc = "Electronic properties";

                    var x = (i, name, desc, values[i]);

                    result.Add(x);
                }

                return result;
            }
        }

        /// <summary>
        /// This function compute the hmoment based on Eisenberg, D., Weiss, R. M., & Terwilliger, T. C. (1984). Hydriphobic moment is a quantitative measure of the amphiphilicity perpendicular to the axis of any periodic peptide structure, such as the a-helix or b-sheet. It can be calculated for an amino acid sequence of N residues and their associated hydrophobicities Hn.
        ///
        /// Eisenberg, D., Weiss, R. M., & Terwilliger, T. C. (1984). The hydrophobic moment detects periodicity in protein hydrophobicity. Proceedings of the National Academy of Sciences, 81(1), 140-144.
        ///
        /// https://rdrr.io/cran/Peptides/man/hmoment.html
        /// </summary>
        /// <param name="seq">An amino-acids sequence</param>
        /// <param name="angle">A protein rotational angle (Suggested: a-helix = 100, b-sheet=160)</param>
        /// <param name="window">A sequence fraction length</param>
        /// <returns>The computed maximal hydrophobic moment (uH) for a given amino-acids sequence</returns>
        public static double hmoment(REngine engine, string seq, int angle = 100, int window = 11)
        {
            lock (engine_lock)
            {
                if (string.IsNullOrWhiteSpace(seq)) return 0;

                var f = nameof(hmoment);
                var k = Key;
                var v = $"{f}_v{k}";

                var evaluate = engine.Evaluate($"{v} <- {f}(seq = \"{seq}\", angle = {angle}, window = {window})");
                var values = engine.Evaluate($"{v}[[1]]").AsNumeric();
                var rm = engine.Evaluate($"rm({v})");

                var result = values.First();


                return result;
            }
        }


        public static string[] hydrophobicity_scales = new string[]
        {
            "Aboderin", "AbrahamLeo", "Argos", "BlackMould", "BullBreese", "Casari", "Chothia", "Cid", "Cowan3.4", "Cowan7.5", "Eisenberg", "Engelman", "Fasman", "Fauchere", "Goldsack", "Guy", "HoppWoods", "Janin", "Jones", "Juretic", "Kidera", "Kuhn", "KyteDoolittle", "Levitt", "Manavalan", "Miyazawa", "Parker", "Ponnuswamy", "Prabhakaran", "Rao", "Rose", "Roseman", "Sweet", "Tanford", "Welling", "Wilson", "Wolfenden", "Zimmerman", "interfaceScale_pH8", "interfaceScale_pH2", "octanolScale_pH8", "octanolScale_pH2", "oiScale_pH8", "oiScale_pH2"
        };

        /// <summary>
        /// This function calculates the GRAVY hydrophobicity index of an amino acids sequence using one of the 38 scales from different sources.
        ///
        /// The hydrophobicity is an important stabilization force in protein folding; this force changes depending on the solvent in which the protein is found. The hydrophobicity index is calculated adding the hydrophobicity of individual amino acids and dividing this value by the length of the sequence.
        ///
        /// Aboderin, A. A. (1971). An empirical hydrophobicity scale for alpha-amino-acids and some of its applications. International Journal of Biochemistry, 2(11), 537-544.
        /// Abraham D.J., Leo A.J. Hydrophobicity (delta G1/2 cal). Proteins: Structure, Function and Genetics 2:130-152(1987).
        /// Argos, P., Rao, J. K., &amp; Hargrave, P. A. (1982). Structural Prediction of Membrane-Bound Proteins. European Journal of Biochemistry, 128(2-3), 565-575.
        /// Black S.D., Mould D.R. Hydrophobicity of physiological L-alpha amino acids. Anal. Biochem. 193:72-82(1991).
        /// Bull H.B., Breese K. Hydrophobicity (free energy of transfer to surface in kcal/mole). Arch. Biochem. Biophys. 161:665-670(1974).
        /// Casari, G., &amp; Sippl, M. J. (1992). Structure-derived hydrophobic potential: hydrophobic potential derived from X-ray structures of globular proteins is able to identify native folds. Journal of molecular biology, 224(3), 725-732.
        /// Chothia, C. (1976). The nature of the accessible and buried surfaces in proteins. Journal of molecular biology, 105(1), 1-12.
        /// Cid, H., Bunster, M., Canales, M., &amp; Gazitua, F. (1992). Hydrophobicity and structural classes in proteins. Protein engineering, 5(5), 373-375.
        /// Cowan R., Whittaker R.G. Hydrophobicity indices at pH 3.4 determined by HPLC. Peptide Research 3:75-80(1990).
        /// Cowan R., Whittaker R.G. Hydrophobicity indices at pH 7.5 determined by HPLC. Peptide Research 3:75-80(1990).
        /// Eisenberg D., Schwarz E., Komarony M., Wall R. Normalized consensus hydrophobicity scale. J. Mol. Biol. 179:125-142(1984).
        /// Engelman, D. M., Steitz, T. A., &amp; Goldman, A. (1986). Identifying nonpolar transbilayer helices in amino acid sequences of membrane proteins. Annual review of biophysics and biophysical chemistry, 15(1), 321-353.
        /// Fasman, G. D. (Ed.). (1989). Prediction of protein structure and the principles of protein conformation. Springer.
        /// Fauchere J.-L., Pliska V.E. Hydrophobicity scale (pi-r). Eur. J. Med. Chem. 18:369-375(1983).
        /// Goldsack, D. E., &amp; Chalifoux, R. C. (1973). Contribution of the free energy of mixing of hydrophobic side chains to the stability of the tertiary structure of proteins. Journal of theoretical biology, 39(3), 645-651.
        /// Guy H.R. Hydrophobicity scale based on free energy of transfer (kcal/mole). Biophys J. 47:61-70(1985).
        /// Hopp T.P., Woods K.R. Hydrophilicity. Proc. Natl. Acad. Sci. U.S.A. 78:3824-3828(1981).
        /// Janin J. Free energy of transfer from inside to outside of a globular protein. Nature 277:491-492(1979).
        /// Jones, D. D. (1975). Amino acid properties and side-chain orientation in proteins: a cross correlation approach. Journal of theoretical biology, 50(1), 167-183.
        /// Juretic, D., Lucic, B., Zucic, D., &amp; Trinajstic, N. (1998). Protein transmembrane structure: recognition and prediction by using hydrophobicity scales through preference functions. Theoretical and computational chemistry, 5, 405-445.
        /// Kidera, A., Konishi, Y., Oka, M., Ooi, T., &amp; Scheraga, H. A. (1985). Statistical analysis of the physical properties of the 20 naturally occurring amino acids. Journal of Protein Chemistry, 4(1), 23-55.
        /// Kuhn, L. A., Swanson, C. A., Pique, M. E., Tainer, J. A., &amp; Getzoff, E. D. (1995). Atomic and residue hydrophilicity in the context of folded protein structures. Proteins: Structure, Function, and Bioinformatics, 23(4), 536-547.
        /// Kyte J., Doolittle R.F. Hydropathicity. J. Mol. Biol. 157:105-132(1982).
        /// Levitt, M. (1976). A simplified representation of protein conformations for rapid simulation of protein folding. Journal of molecular biology, 104(1), 59-107.
        /// Manavalan P., Ponnuswamy Average surrounding hydrophobicity. P.K. Nature 275:673-674(1978).
        /// Miyazawa S., Jernigen R.L. Hydrophobicity scale (contact energy derived from 3D data). Macromolecules 18:534-552(1985).
        /// Parker J.M.R., Guo D., Hodges R.S. Hydrophilicity scale derived from HPLC peptide retention times. Biochemistry 25:5425-5431(1986).
        /// Ponnuswamy, P. K. (1993). Hydrophobic charactesristics of folded proteins. Progress in biophysics and molecular biology, 59(1), 57-103.
        /// Prabhakaran, M. (1990). The distribution of physical, chemical and conformational properties in signal and nascent peptides. Biochem. J, 269, 691-696.
        /// Rao M.J.K., Argos P. Membrane buried helix parameter. Biochim. Biophys. Acta 869:197-214(1986).
        /// Rose G.D., Geselowitz A.R., Lesser G.J., Lee R.H., Zehfus M.H. Mean fractional area loss (f) [average area buried/standard state area]. Science 229:834-838(1985)
        /// Roseman M.A. Hydrophobicity scale (pi-r). J. Mol. Biol. 200:513-522(1988).
        /// Sweet R.M., Eisenberg D. Optimized matching hydrophobicity (OMH). J. Mol. Biol. 171:479-488(1983).
        /// Tanford C. Hydrophobicity scale (Contribution of hydrophobic interactions to the stability of the globular conformation of proteins). J. Am. Chem. Soc. 84:4240-4274(1962).
        /// Welling G.W., Weijer W.J., Van der Zee R., Welling-Wester S. Antigenicity value X 10. FEBS Lett. 188:215-218(1985).
        /// Wilson K.J., Honegger A., Stotzel R.P., Hughes G.J. Hydrophobic constants derived from HPLC peptide retention times. Biochem. J. 199:31-41(1981).
        /// Wolfenden R.V., Andersson L., Cullis P.M., Southgate C.C.F. Hydration potential (kcal/mole) at 25C. Biochemistry 20:849-855(1981).
        /// Zimmerman, J. M., Eliezer, N., &amp; Simha, R. (1968). The characterization of amino acid sequences in proteins by statistical methods. Journal of theoretical biology, 21(2), 170-201.
        /// Nakai, K., Kidera, A., and Kanehisa, M.; Cluster analysis of amino acid indices for prediction of protein structure and function. Protein Eng. 2, 93-100 (1988).
        /// Tomii, K. and Kanehisa, M.; Analysis of amino acid indices and mutation matrices for sequence comparison and structure prediction of proteins. Protein Eng. 9, 27-36 (1996).
        /// Kawashima, S., Ogata, H., and Kanehisa, M.; AAindex: amino acid index database. Nucleic Acids Res. 27, 368-369 (1999).
        /// Kawashima, S. and Kanehisa, M.; AAindex: amino acid index database. Nucleic Acids Res. 28, 374 (2000).
        /// Kawashima, S., Pokarowski, P., Pokarowska, M., Kolinski, A., Katayama, T., and Kanehisa, M.; AAindex: amino acid index database, progress report 2008. Nucleic Acids Res. 36, D202-D205 (2008).
        /// White, Stephen (2006-06-29). "Experimentally Determined Hydrophobicity Scales". University of California, Irvine. Retrieved 2017-05-25
        ///
        /// https://rdrr.io/cran/Peptides/man/hydrophobicity.html
        /// </summary>
        /// <param name="seq">An amino-acids sequence</param>
        /// <param name="scale">A character string specifying the hydophobicity scale to be used; must be one of "Aboderin", "AbrahamLeo", "Argos", "BlackMould", "BullBreese", "Casari", "Chothia", "Cid", "Cowan3.4", "Cowan7.5", "Eisenberg", "Engelman", "Fasman", "Fauchere", "Goldsack", "Guy", "HoppWoods", "Janin", "Jones", "Juretic", "Kidera", "Kuhn", "KyteDoolittle", "Levitt", "Manavalan", "Miyazawa", "Parker", "Ponnuswamy", "Prabhakaran", "Rao", "Rose", "Roseman", "Sweet", "Tanford", "Welling", "Wilson", "Wolfenden", "Zimmerman", "interfaceScale_pH8", "interfaceScale_pH2", "octanolScale_pH8", "octanolScale_pH2", "oiScale_pH8" or "oiScale_pH2".</param>
        /// <returns>The computed GRAVY index for a given amino-acid sequence</returns>
        public static double hydrophobicity(REngine engine, string seq, string scale = "KyteDoolittle")
        {
            lock (engine_lock)
            {
                if (!hydrophobicity_scales.Contains(scale)) throw new Exception();

                var f = nameof(hydrophobicity);
                var k = Key;
                var v = $"{f}_v{k}";

                var eval_cmd = $"{v} <- {f}({nameof(seq)} = \"{seq}\", {nameof(scale)} = \"{scale}\")";
                var evaluate = engine.Evaluate(eval_cmd);
                var values = engine.Evaluate($"{v}[[1]]").AsNumeric();
                var rm = engine.Evaluate($"rm({v})");

                var result = values.First();


                return result;
            }
        }


        /// <summary>
        /// This function calculates the instability index proposed by Guruprasad (1990). This index predicts the stability of a protein based on its amino acid composition, a protein whose instability index is smaller than 40 is predicted as stable, a value above 40 predicts that the protein may be unstable.
        ///
        /// Guruprasad K, Reddy BV, Pandit MW (1990). "Correlation between stability of a protein and its dipeptide composition: a novel approach for predicting in vivo stability of a protein from its primary sequence". Protein Eng. 4 (2): 155 - 61. doi:10.1093/protein/4.2.155
        ///
        /// https://rdrr.io/cran/Peptides/man/instaIndex.html
        /// </summary>
        /// <param name="seq">An amino-acids sequence</param>
        /// <returns>The computed instability index for a given amino-acids sequence</returns>
        public static double instaIndex(REngine engine, string seq)
        {
            lock (engine_lock)
            {
                var f = nameof(instaIndex);
                var k = Key;
                var v = $"{f}_v{k}";

                if (seq.Length <= 2) return 0;

                var eval_cmd = $"{v} <- {f}({nameof(seq)} = \"{seq}\")";
                var evaluate = engine.Evaluate(eval_cmd);
                var values = engine.Evaluate($"{v}[[1]]").AsNumeric();
                var rm = engine.Evaluate($"rm({v})");

                var result = values.First();


                return result;
            }
        }

        /// <summary>
        /// The Kidera Factors were originally derived by applying multivariate analysis to 188 physical properties of the 20 amino acids and using dimension reduction techniques. This function calculates the average of the ten Kidera factors for a protein sequence.
        ///
        /// Kidera, A., Konishi, Y., Oka, M., Ooi, T., & Scheraga, H. A. (1985). Statistical analysis of the physical properties of the 20 naturally occurring amino acids. Journal of Protein Chemistry, 4(1), 23-55.
        ///
        /// https://rdrr.io/cran/Peptides/man/kideraFactors.html
        /// </summary>
        /// <param name="seq">An amino-acids sequence</param>
        /// <returns>
        /// A list with the average of the ten Kidera factors. The first four factors are essentially pure physical properties; the remaining six factors are superpositions of several physical properties, and are labelled for convenience by the name of the most heavily weighted component.
        /// 
        /// KF1: Helix/bend preference,
        /// KF2: Side-chain size,
        /// KF3: Extended structure preference,
        /// KF4: Hydrophobicity,
        /// KF5: Double-bend preference,
        /// KF6: Partial specific volume,
        /// KF7: Flat extended preference,
        /// KF8: Occurrence in alpha region,
        /// KF9: pK-C,
        /// KF10: Surrounding hydrophobicity
        /// </returns>
        public static List<(int index, string name, string description, double value)> kideraFactors(REngine engine, string seq)
        {
            lock (engine_lock)
            {
                var f = nameof(kideraFactors);
                var k = Key;
                var v = $"{f}_v{k}";
                var ai = "[[1]]";

                var evaluate = engine.Evaluate($"{v} <- {f}(seq = \"{seq}\")");
                var names = engine.Evaluate($"names({v}{ai})").AsCharacter();
                var dimnames = engine.Evaluate($"dimnames({v}{ai})");
                var rownames = engine.Evaluate($"rownames({v}{ai})");
                var colnames = engine.Evaluate($"colnames({v}{ai})");
                var values = engine.Evaluate($"{v}{ai}").AsNumeric();
                var rm = engine.Evaluate($"rm({v})");

                var result = new List<(int index, string name, string description, double value)>();

                for (var i = 0; i < names.Length; i++)
                {
                    var name = names[i];
                    var desc = "";
                    if (name == "KF1") desc = "Helix/bend preference";
                    else if (name == "KF2") desc = "Side-chain size";
                    else if (name == "KF3") desc = "Extended structure preference";
                    else if (name == "KF4") desc = "Hydrophobicity";
                    else if (name == "KF5") desc = "Double-bend preference";
                    else if (name == "KF6") desc = "Partial specific volume";
                    else if (name == "KF7") desc = "Flat extended preference";
                    else if (name == "KF8") desc = "Occurrence in alpha region";
                    else if (name == "KF9") desc = "pK-C";
                    else if (name == "KF10") desc = "Surrounding hydrophobicity";

                    var x = (i, name, desc, values[i]);

                    result.Add(x);
                }

                return result;
            }
        }

        /// <summary>
        /// This function counts the number of amino acids in a protein sequence
        ///
        /// All proteins are formed by linear chains of small residues known as amino acids attached to each other by peptide bonds. The function lengthpep counts the number of amino acids in a sequence and returns a vector with the count for each peptide used as argument.
        ///
        /// https://rdrr.io/cran/Peptides/man/lengthpep.html
        /// </summary>
        /// <param name="seq">An amino-acids sequence</param>
        /// <returns></returns>
        public static double lengthpep(REngine engine, string seq)
        {
            lock (engine_lock)
            {
                var f = nameof(lengthpep);
                var k = Key;
                var v = $"{f}_v{k}";

                var eval_cmd = $"{v} <- {f}({nameof(seq)} = \"{seq}\")";
                var evaluate = engine.Evaluate(eval_cmd);
                var values = engine.Evaluate($"{v}[[1]]").AsNumeric();
                var rm = engine.Evaluate($"rm({v})");

                var result = values.First();


                return result;
            }
        }


        /// <summary>
        /// This function calculates the theoretical class of a protein sequence based on the relationship between the hydrophobic moment and hydrophobicity scale proposed by Eisenberg (1984).
        ///
        /// Eisenberg et al. (1982) found a correlation between hydrophobicity and hydrophobic moment that defines the protein section as globular, transmembrane or superficial. The function calculates the hydrophobicity (H) and hydrophobic moment (uH) based on the standardized scale of Eisenberg (1984) using windows of 11 amino acids for calculate the theoretical fragment type.
        /// 
        /// Eisenberg, David. "Three-dimensional structure of membrane and surface proteins." Annual review of biochemistry 53.1 (1984): 595-623.
        /// 
        /// D. Eisenberg, R. M. Weiss, and T. C. Terwilliger. The helical hydrophobic moment: A measure of the amphiphilicity of a helix. Nature, 299(5881):371-374, 1982. [p7, 8]
        ///
        /// https://rdrr.io/cran/Peptides/man/membpos.html
        /// </summary>
        /// <param name="seq">An amino-acids sequence</param>
        /// <param name="angle">A protein rotational angle</param>
        /// <returns>A data frame for each sequence given with the calculated class for each window of eleven amino-acids</returns>
        public static List<(int row, string Pwp, double H, double uH, double MembPosValue, string MembPos)> membpos(REngine engine, string seq, int angle = 100)
        {
            lock (engine_lock)
            {
                if (string.IsNullOrWhiteSpace(seq)) return new List<(int row, string Pwp, double H, double uH, double MembPosValue, string MembPos)>();

                var f = nameof(membpos);
                var k = Key;
                var v = $"{f}_v{k}";
                var ai = "[[1]]";

                var eval_cmd = $"{v} <- {f}({nameof(seq)} = \"{seq}\", {nameof(angle)} = {angle})";

                var evaluate = engine.Evaluate(eval_cmd);

                var names = engine.Evaluate($"names({v}{ai})").AsCharacter();
                var dimnames = engine.Evaluate($"dimnames({v}{ai})").AsCharacter();
                var rownames = engine.Evaluate($"rownames({v}{ai})").AsCharacter();
                var colnames = engine.Evaluate($"colnames({v}{ai})").AsCharacter();

                var values_Pep = engine.Evaluate($"{v}{ai}[[\"Pep\"]]").AsCharacter();
                var values_H = engine.Evaluate($"{v}{ai}[[\"H\"]]").AsNumeric();
                var values_uH = engine.Evaluate($"{v}{ai}[[\"uH\"]]").AsNumeric();
                var values_MembPos = engine.Evaluate($"{v}{ai}[[\"MembPos\"]]").AsCharacter();
                var rm = engine.Evaluate($"rm({v})");

                var list = new List<(int row, string Pwp, double H, double uH, double MembPosValue, string MembPos)>();

                for (var r = 0; r < rownames.Length; r++)
                {
                    var MembPos = values_MembPos[r];
                    var MembPosValue = 0.5;

                    if (MembPos == "Globular") MembPosValue = 0.5;
                    else if (MembPos == "Surface") MembPosValue = 1.0;
                    else if (MembPos == "Transmembrane") MembPosValue = 0.0;


                    var x = (r, values_Pep[r], values_H[r], values_uH[r], MembPosValue, MembPos);

                    list.Add(x);
                }

                return list;
            }
        }

        /// <summary>
        /// MS-WHIM scores were derived from 36 electrostatic potential properties derived from the three-dimensional structure of the 20 natural amino acids
        ///
        /// Zaliani, A., & Gancia, E. (1999). MS-WHIM scores for amino acids: a new 3D-description for peptide QSAR and QSPR studies. Journal of chemical information and computer sciences, 39(3), 525-533.
        ///
        /// https://rdrr.io/cran/Peptides/man/mswhimScores.html
        /// </summary>
        /// <param name="seq">An amino-acids sequence</param>
        /// <returns>The computed average of MS-WHIM scores of all the amino acids in the corresponding peptide sequence.</returns>
        public static List<(int index, string name, double value)> mswhimScores(REngine engine, string seq)
        {
            lock (engine_lock)
            {
                var f = nameof(mswhimScores);
                var k = Key;
                var v = $"{f}_v{k}";
                var ai = "[[1]]";

                var evaluate = engine.Evaluate($"{v} <- {f}(seq = \"{seq}\")");
                var names = engine.Evaluate($"names({v}{ai})").AsCharacter();
                var dimnames = engine.Evaluate($"dimnames({v}{ai})");
                var rownames = engine.Evaluate($"rownames({v}{ai})");
                var colnames = engine.Evaluate($"colnames({v}{ai})");
                var values = engine.Evaluate($"{v}{ai}").AsNumeric();
                var rm = engine.Evaluate($"rm({v})");

                var result = new List<(int index, string name, double value)>();

                for (var i = 0; i < names.Length; i++)
                {
                    var name = names[i];


                    var x = (i, name, values[i]);

                    result.Add(x);
                }

                return result;
            }
        }

        /// <summary>
        /// This function calculates the molecular weight of a protein sequence. It is calculated as the sum of the mass of each amino acid using the scale available on Compute pI/Mw tool.
        ///
        /// The molecular weight is the sum of the masses of each atom constituting a molecule. The molecular weight is directly related to the length of the amino acid sequence and is expressed in units called daltons (Da). In Peptides the function mw computes the molecular weight using the same formulas and weights as ExPASy's "compute pI/mw" tool (Gasteiger et al., 2005).
        ///
        /// The formula and amino acid scale are the same available on ExPASy Compute pI/Mw tool: http://web.expasy.org/compute_pi/
        ///
        /// Gasteiger, E., Hoogland, C., Gattiker, A., Wilkins, M. R., Appel, R. D., & Bairoch, A. (2005). Protein identification and analysis tools on the ExPASy server. In The proteomics protocols handbook (pp. 571-607). Humana Press. Chicago
        ///
        /// https://rdrr.io/cran/Peptides/man/mw.html
        /// </summary>
        /// <param name="seq"></param>
        /// <returns></returns>
        public static double mw(REngine engine, string seq, bool monoisotopic = false)
        {
            lock (engine_lock)
            {


                var f = nameof(mw);
                var k = Key;
                var v = $"{f}_v{k}";

                var eval_cmd = $"{v} <- {f}({nameof(seq)} = \"{seq}\", {nameof(monoisotopic)} = {monoisotopic.ToString().ToUpperInvariant()})";
                var evaluate = engine.Evaluate(eval_cmd);
                var values = engine.Evaluate($"{v}[[1]]").AsNumeric();
                var rm = engine.Evaluate($"rm({v})");

                var result = values.First();



                return result;
            }
        }


        /// <summary>
        /// The isoelectric point (pI), is the pH at which a particular molecule or surface carries no net electrical charge.
        ///
        /// The isoelectric point (pI) is the pH at which the net charge of the protein is equal to 0. It is a variable that affects the solubility of the peptides under certain conditions of pH. When the pH of the solvent is equal to the pI of the protein, it tends to precipitate and loose its biological function.
        ///
        /// https://rdrr.io/cran/Peptides/man/pI.html
        /// </summary>
        /// <param name="seq">An amino-acids sequence</param>
        /// <param name="pKscale">A character string specifying the pK scale to be used; must be one of "Bjellqvist", "EMBOSS", "Murray", "Sillero", "Solomon", "Stryer", "Lehninger", "Dawson" or "Rodwell"</param>
        /// <returns></returns>
        public static double pI(REngine engine, string seq, string pKscale = "EMBOSS")
        {
            lock (engine_lock)
            {
                if (!pKscales.Contains(pKscale)) throw new Exception();

                var f = nameof(pI);
                var k = Key;
                var v = $"{f}_v{k}";

                var eval_cmd = $"{v} <- {f}({nameof(seq)} = \"{seq}\", {nameof(pKscale)} = \"{pKscale}\")";
                var evaluate = engine.Evaluate(eval_cmd);
                var values = engine.Evaluate($"{v}[[1]]").AsNumeric();
                var rm = engine.Evaluate($"rm({v})");

                var result = values.First();


                return result;
            }
        }

        /// <summary>
        /// The ProtFP descriptor set was constructed from a large initial selection of indices obtained from the AAindex database for all 20 naturally occurring amino acids.
        ///
        /// van Westen, G. J., Swier, R. F., Wegner, J. K., IJzerman, A. P., van Vlijmen, H. W., & Bender, A. (2013). Benchmarking of protein descriptor sets in proteochemometric modeling (part 1): comparative study of 13 amino acid descriptor sets. Journal of cheminformatics, 5(1), 41.
        ///
        /// https://rdrr.io/cran/Peptides/man/protFP.html
        /// </summary>
        /// <param name="seq">An amino-acids sequence</param>
        /// <returns>The computed average of protFP descriptors of all the amino acids in the corresponding peptide sequence.</returns>
        public static List<(int index, string name, double value)> protFP(REngine engine, string seq)
        {
            lock (engine_lock)
            {
                var f = nameof(protFP);
                var k = Key;
                var v = $"{f}_v{k}";
                var ai = "[[1]]";

                var evaluate = engine.Evaluate($"{v} <- {f}(seq = \"{seq}\")");
                var names = engine.Evaluate($"names({v}{ai})").AsCharacter();
                var dimnames = engine.Evaluate($"dimnames({v}{ai})");
                var rownames = engine.Evaluate($"rownames({v}{ai})");
                var colnames = engine.Evaluate($"colnames({v}{ai})");
                var values = engine.Evaluate($"{v}{ai}").AsNumeric();
                var rm = engine.Evaluate($"rm({v})");

                var result = new List<(int index, string name, double value)>();

                for (var i = 0; i < names.Length; i++)
                {
                    var name = names[i];


                    var x = (i, name, values[i]);

                    result.Add(x);
                }

                return result;
            }
        }

        /// <summary>
        /// ST-scales were proposed by Yang et al, taking 827 properties into account which are mainly constitutional, topological, geometrical, hydrophobic, elec- tronic, and steric properties of a total set of 167 AAs.
        ///
        /// Yang, L., Shu, M., Ma, K., Mei, H., Jiang, Y., & Li, Z. (2010). ST-scale as a novel amino acid descriptor and its application in QSAM of peptides and analogues. Amino acids, 38(3), 805-816.
        ///
        /// https://rdrr.io/cran/Peptides/man/stScales.html
        /// </summary>
        /// <param name="seq">An amino-acids sequence</param>
        /// <returns>The computed average of ST-scales of all the amino acids in the corresponding peptide sequence.</returns>
        public static List<(int index, string name, double value)> stScales(REngine engine, string seq)
        {
            lock (engine_lock)
            {
                var f = nameof(stScales);
                var k = Key;
                var v = $"{f}_v{k}";
                var ai = "[[1]]";

                var evaluate = engine.Evaluate($"{v} <- {f}(seq = \"{seq}\")");
                var names = engine.Evaluate($"names({v}{ai})").AsCharacter();
                var dimnames = engine.Evaluate($"dimnames({v}{ai})");
                var rownames = engine.Evaluate($"rownames({v}{ai})");
                var colnames = engine.Evaluate($"colnames({v}{ai})");
                var values = engine.Evaluate($"{v}{ai}").AsNumeric();
                var rm = engine.Evaluate($"rm({v})");

                var result = new List<(int index, string name, double value)>();

                for (var i = 0; i < names.Length; i++)
                {
                    var name = names[i];


                    var x = (i, name, values[i]);

                    result.Add(x);
                }

                return result;
            }
        }

        /// <summary>
        /// T-scales are based on 67 common topological descriptors of 135 amino acids. These topological descriptors are based on the connectivity table of amino acids alone, and to not explicitly consider 3D properties of each structure.
        ///
        /// Tian F, Zhou P, Li Z: T-scale as a novel vector of topological descriptors for amino acids and its application in QSARs of peptides. J Mol Struct. 2007, 830: 106-115. 10.1016/j.molstruc.2006.07.004.
        ///
        /// https://rdrr.io/cran/Peptides/man/tScales.html
        /// </summary>
        /// <param name="seq">An amino-acids sequence</param>
        /// <returns>The computed average of T-scales of all the amino acids in the corresponding peptide sequence.</returns>
        public static List<(int index, string name, double value)> tScales(REngine engine, string seq)
        {
            lock (engine_lock)
            {
                var f = nameof(tScales);
                var k = Key;
                var v = $"{f}_v{k}";
                var ai = "[[1]]";

                var evaluate = engine.Evaluate($"{v} <- {f}(seq = \"{seq}\")");
                var names = engine.Evaluate($"names({v}{ai})").AsCharacter();
                var dimnames = engine.Evaluate($"dimnames({v}{ai})");
                var rownames = engine.Evaluate($"rownames({v}{ai})");
                var colnames = engine.Evaluate($"colnames({v}{ai})");
                var values = engine.Evaluate($"{v}{ai}").AsNumeric();
                var rm = engine.Evaluate($"rm({v})");

                var result = new List<(int index, string name, double value)>();

                for (var i = 0; i < names.Length; i++)
                {
                    var name = names[i];


                    var x = (i, name, values[i]);

                    result.Add(x);
                }

                return result;
            }
        }

        /// <summary>
        /// VHSE-scales (principal components score Vectors of Hydrophobic, Steric, and Electronic properties), is derived from principal components analysis (PCA) on independent families of 18 hydrophobic properties, 17 steric properties, and 15 electronic properties, respectively, which are included in total 50 physicochemical variables of 20 coded amino acids.
        ///
        /// Mei, H. U., Liao, Z. H., Zhou, Y., & Li, S. Z. (2005). A new set of amino acid descriptors and its application in peptide QSARs. Peptide Science, 80(6), 775-786.
        ///
        /// https://rdrr.io/cran/Peptides/man/vhseScales.html
        /// </summary>
        /// <param name="seq">An amino-acids sequence</param>
        /// <returns>
        /// The computed average of VHSE-scales of all the amino acids in the corresponding peptide sequence. Each VSHE-scale represent an amino-acid property as follows:
        /// VHSE1 and VHSE2: Hydrophobic properties
        /// VHSE3 and VHSE4: Steric properties
        /// VHSE5 to VHSE8: Electronic properties
        /// </returns>
        public static List<(int index, string name, string description, double value)> vhseScales(REngine engine, string seq)
        {
            lock (engine_lock)
            {
                var f = nameof(vhseScales);
                var k = Key;
                var v = $"{f}_v{k}";
                var ai = "[[1]]";

                var evaluate = engine.Evaluate($"{v} <- {f}(seq = \"{seq}\")");
                var names = engine.Evaluate($"names({v}{ai})").AsCharacter();
                var dimnames = engine.Evaluate($"dimnames({v}{ai})");
                var rownames = engine.Evaluate($"rownames({v}{ai})");
                var colnames = engine.Evaluate($"colnames({v}{ai})");
                var values = engine.Evaluate($"{v}{ai}").AsNumeric();
                var rm = engine.Evaluate($"rm({v})");

                var result = new List<(int index, string name, string description, double value)>();

                for (var i = 0; i < names.Length; i++)
                {
                    var name = names[i];

                    var desc = "";
                    if (name == "VHSE1") desc = "Hydrophobic properties";
                    else if (name == "VHSE2") desc = "Hydrophobic properties";
                    else if (name == "VHSE3") desc = "Steric properties";
                    else if (name == "VHSE4") desc = "Steric properties";
                    else if (name == "VHSE5") desc = "Electronic properties";
                    else if (name == "VHSE6") desc = "Electronic properties";
                    else if (name == "VHSE7") desc = "Electronic properties";
                    else if (name == "VHSE8") desc = "Electronic properties";

                    var x = (i, name, desc, values[i]);

                    result.Add(x);
                }

                return result;
            }
        }


        /// <summary>
        /// Z-scales are based on physicochemical properties of the AAs including NMR data and thin-layer chromatography (TLC) data.
        ///
        /// Sandberg M, Eriksson L, Jonsson J, Sjostrom M, Wold S: New chemical descriptors relevant for the design of biologically active peptides. A multivariate characterization of 87 amino acids. J Med Chem 1998, 41:2481-2491.
        ///
        /// https://rdrr.io/cran/Peptides/man/zScales.html
        /// </summary>
        /// <param name="seq">An amino-acids sequence</param>
        /// <returns>
        /// The computed average of Z-scales of all the amino acids in the corresponding peptide sequence. Each Z scale represent an amino-acid property as follows:
        /// Z1: Lipophilicity
        /// Z2: Steric properties (Steric bulk/Polarizability)
        /// Z3: Electronic properties (Polarity / Charge)
        /// Z4 and Z5: They relate electronegativity, heat of formation, electrophilicity and hardness.
        /// </returns>
        public static List<(int index, string name, string description, double value)> zScales(REngine engine, string seq)
        {
            lock (engine_lock)
            {
                var f = nameof(zScales);
                var k = Key;
                var v = $"{f}_v{k}";
                var ai = "[[1]]";

                var evaluate = engine.Evaluate($"{v} <- {f}(seq = \"{seq}\")");
                var names = engine.Evaluate($"names({v}{ai})").AsCharacter();
                var dimnames = engine.Evaluate($"dimnames({v}{ai})");
                var rownames = engine.Evaluate($"rownames({v}{ai})");
                var colnames = engine.Evaluate($"colnames({v}{ai})");
                var values = engine.Evaluate($"{v}{ai}").AsNumeric();
                var rm = engine.Evaluate($"rm({v})");

                var result = new List<(int index, string name, string description, double value)>();

                for (var i = 0; i < names.Length; i++)
                {
                    var name = names[i];
                    var desc = "";

                    if (name == "Z1") desc = "Lipophilicity";
                    else if (name == "Z2") desc = "Steric properties (Steric bulk/Polarizability)";
                    else if (name == "Z3") desc = "Electronic properties (Polarity / Charge)";
                    else if (name == "Z4") desc = "Electronegativity, heat of formation, electrophilicity and hardness";
                    else if (name == "Z5") desc = "Electronegativity, heat of formation, electrophilicity and hardness";

                    var x = (i, name, desc, values[i]);

                    result.Add(x);
                }

                return result;
            }
        }


    }
}
