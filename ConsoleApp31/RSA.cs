using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;

namespace Coil_DHC_Dataset_Maker
{
    public class RSA
    {
        public string pdb_id;
        public char algo;

        public char amino_acid;
        public char chain_id;
        public int res_num;

        public int array_index;

        public double all_atoms_abs;
        public double all_atoms_rel; 

        public double total_side_abs;
        public double total_side_rel;

        public double main_chain_abs;
        public double main_chain_rel;


        public double non_polar_abs;
        public double non_polar_rel;

        public double all_polar_abs;
        public double all_polar_rel;

        



        public static List<RSA> load(List<string> files)
        {
            var rsa_list = new List<RSA>();

            foreach (var file in files)
            {
                var pdb_id = Path.GetFileNameWithoutExtension(file).Substring(0, 4).ToUpperInvariant();

                var algo = Path.GetExtension(file)[1];

                var lines = File.ReadAllLines(file).ToList();

                lines = lines.Where(a => a.StartsWith("RES ")).ToList();


                var array_index = 0;
                foreach (var line in lines)
                {
                    var split = line.Split(new char[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries).ToList();

                    if (split[2].Any(char.IsDigit) && split[2].Any(char.IsLetter)) { var a = split[2][0]; var b = split[2].Substring(1); split.RemoveAt(2); split.Insert(2, "" + a); split.Insert(3, b); }

                    var i = 0;//skip first
                    var rsa = new RSA()
                    {
                        pdb_id = pdb_id,
                        algo = algo,
                        array_index = ++array_index,
                        amino_acid = Atom.Aa3To1(split[++i]),

                        chain_id = split[++i][0],
                        res_num = int.Parse(split[++i]),

                        all_atoms_abs = split[++i].try_parse_double(),
                        all_atoms_rel = split[++i].try_parse_double(),
                        total_side_abs =split[++i].try_parse_double(),
                        total_side_rel =split[++i].try_parse_double(),
                        main_chain_abs =split[++i].try_parse_double(),
                        main_chain_rel =split[++i].try_parse_double(),
                        non_polar_abs = split[++i].try_parse_double(),
                        non_polar_rel = split[++i].try_parse_double(),
                        all_polar_abs = split[++i].try_parse_double(),
                        all_polar_rel = split[++i].try_parse_double(),
                    };

                    rsa_list.Add(rsa);
                }
            }

            return rsa_list;
        }
    }
}
