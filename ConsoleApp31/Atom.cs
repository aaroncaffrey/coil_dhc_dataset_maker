using System;
using System.Collections.Generic;
using System.Data.Odbc;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Linq.Expressions;
using System.Net;
using System.Net.Configuration;
using System.Reflection;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace Coil_DHC_Dataset_Maker
{
    
    public class Atom
    {
        public string uniprot_sequence;

        public string pdb_id;
        public int model_index;
        public char chain_id;
        public string atom_type;
        public char element;
        public char amino_acid;
        public int residue_index = -1;
        public int serial_index = -1;

        public int array_index = -1;

        public int master_index = -1;

        public char i_code;

        public char multimer_dssp = ' ';
        public char multimer_dssp3 = ' ';

        public char monomer_dssp = ' ';
        public char monomer_dssp3 = ' ';

        public char monomer_stride = ' ';
        public char monomer_stride3 = ' ';

        public char multimer_stride = ' ';
        public char multimer_stride3 = ' ';


        public int is_strand_interface_atom = 0;
        public int is_standard_strand_interface_atom = 0;
        public int is_dimorphic_strand_interface_atom = 0;
        public int is_hybrid_strand_interface_atom = 0;

        public char parallelism = '_';
        public char symmetry_mode = '_';
        public char strand_interface_type = '_';

        public double X;
        public double Y;
        public double Z;

        public RSA RSA_L;
        public RSA RSA_S;

        public sable.sable_item sable_item;

        public List<(string database, List<pssm.pssm_entry> pssm_entries)> amino_acid_pssm_unnormalised;
        public List<(string database, List<pssm.pssm_entry> pssm_entries)> amino_acid_pssm_normalised;

        public List<(string format, mpsa_reader.mpsa_line_entry mpsa_entry)> mpsa_entries;

        public (int index, double short_type_score, double long_type_score, double glob_type_score, double anchor2_score) iup_entry;

        public List<Atom> amino_acid_atoms;
        public List<Atom> chain_atoms;

        public Atom amino_acid_master_atom;
        public List<(Atom atom, double distance)> contact_map_intermolecular;
        public List<(Atom atom, double distance)> contact_map_intramolecular;

        public List<(int array_index, char chain_id)> chain_order_in_pdb_file;

        public double foldx_monomer_ala_scan_ddg;

        public List<ring.ring_edge> monomer_ring_edges;
        public List<ring.ring_node> monomer_ring_nodes;


        public bool amino_acid_exists_in_repaired_structure = false;

        public double chain_dna_binding_prob_nr = 0d;
        public double chain_dna_binding_prob_swissprot = 0d;
        public double chain_dna_binding_prob_uniref90 = 0d;

        //public List<(bool repaired, char mutant_res, double ddg)> foldx_position_scan_ddg;

        //META30P




        public static char Aa3To1(string aa)
        {
            var aa_codes = new (string three_letter_code, char one_letter_code)[]
            {
                ("Arg", 'R'), ("Lys", 'K'), ("Asp", 'D'), ("Glu", 'E'), ("Gln", 'Q'), ("Asn", 'N'), ("His", 'H'),
                ("Ser", 'S'), ("Thr", 'T'), ("Tyr", 'Y'), ("Cys", 'C'), ("Met", 'M'), ("Trp", 'W'), ("Ala", 'A'),
                ("Ile", 'I'), ("Leu", 'L'), ("Phe", 'F'), ("Val", 'V'), ("Pro", 'P'), ("Gly", 'G'), ("XAA", 'X'),
                ("UNK", 'X'), ("MSE", 'M')
            };

            var r = aa_codes.FirstOrDefault(a => String.Equals(a.three_letter_code, aa.Trim(), StringComparison.InvariantCultureIgnoreCase)).one_letter_code;

            if (r == default(char))
            {
            }

            return r;
        }

        public static string Aa1To3(char aa)
        {
            if (!char.IsLetter(aa)) return "" + aa;

            var aa_codes = new (string three_letter_code, char one_letter_code)[]
            {
                ("Arg", 'R'), ("Lys", 'K'), ("Asp", 'D'), ("Glu", 'E'), ("Gln", 'Q'), ("Asn", 'N'), ("His", 'H'),
                ("Ser", 'S'), ("Thr", 'T'), ("Tyr", 'Y'), ("Cys", 'C'), ("Met", 'M'), ("Trp", 'W'), ("Ala", 'A'),
                ("Ile", 'I'), ("Leu", 'L'), ("Phe", 'F'), ("Val", 'V'), ("Pro", 'P'), ("Gly", 'G'), ("XAA", 'X'),
                ("UNK", 'X')
            };

            var r = aa_codes.FirstOrDefault(a => a.one_letter_code == aa).three_letter_code;

            if (r == default(string))
            {
                r = "" + aa;
            }

            return r.ToUpperInvariant();
        }

        public static string secondary_structure_state_reduction(string ss)
        {
            if (string.IsNullOrEmpty(ss)) return ss;

            return string.Concat(ss.Select(secondary_structure_state_reduction).ToList());
        }

        public static char secondary_structure_state_reduction(char ss)
        {
            if (ss == 'E' || ss == 'B' || ss == 'e' || ss == 'b') return 'E'; // change E/e, B/b to E
            if (ss == 'G' || ss == 'g' || ss == 'H' || ss == 'h') return 'H'; // change G/g, H/h to H
            //if (ss == 'T' || ss == 'S' || ss == ' ' || ss == '_') return 'C';
            return 'C'; // all others to C
        }

        public static string[] extract_split_pdb_chains(string pdb_id, char? chain_id)
        {
            var use_ram_disk = false;

            var drive_letter = use_ram_disk ? "r" : "c";

            var pdb_in_folder = $@"{drive_letter}:\betastrands_dataset\pdb\";
            var pdb_out_folder = $@"{drive_letter}:\betastrands_dataset\pdb_split\";

            pdb_id = Path.GetFileNameWithoutExtension(pdb_id);
            pdb_id = pdb_id.ToUpperInvariant();

            var pdb_lines = File.ReadAllLines($@"{pdb_in_folder}{pdb_id}.pdb").ToList();
            var endmdl_index = pdb_lines.FindIndex(a => a.StartsWith("ENDMDL"));
            if (endmdl_index > -1) pdb_lines = pdb_lines.Take(endmdl_index).ToList();
            pdb_lines = pdb_lines.Where(a => a.StartsWith("ATOM ")).ToList();

            var chains = pdb_lines.GroupBy(a => a[21]).Where(a => chain_id == null || a.Key == chain_id).ToList();
            chains.ForEach(a =>
            {
                var chain_pdb_file = $@"{pdb_out_folder}{pdb_id}{a.Key}.pdb";
                if (!File.Exists(chain_pdb_file) || new FileInfo(chain_pdb_file).Length == 0)
                {
                    Directory.CreateDirectory(Path.GetDirectoryName(chain_pdb_file));

                    File.WriteAllLines(chain_pdb_file, a.ToList());
                }
            });
            return chains.Select(a => $@"{pdb_out_folder}{pdb_id}{a.Key}.pdb").ToArray();
        }

        public static string extract_split_pdb_chains_res_ids(string pdb_id, char chain_id, int first_res_id, int last_res_id)
        {
            if (first_res_id > last_res_id)
            {
                var first = last_res_id;
                var last = first_res_id;
                first_res_id = first;
                last_res_id = last;
            }

            var use_ram_disk = false;

            var drive_letter = use_ram_disk ? "r" : "c";

            var pdb_in_folder = $@"{drive_letter}:\betastrands_dataset\pdb\";
            var pdb_out_folder = $@"{drive_letter}:\betastrands_dataset\pdb_split\";

            pdb_id = pdb_id.ToUpperInvariant();
            var pdb_lines = File.ReadAllLines($@"{pdb_in_folder}{pdb_id}.pdb").ToList();
            var endmdl_index = pdb_lines.FindIndex(a => a.StartsWith("ENDMDL"));
            if (endmdl_index > -1) pdb_lines = pdb_lines.Take(endmdl_index).ToList();
            pdb_lines = pdb_lines.Where(a => a.StartsWith("ATOM ") && a[21] == chain_id && int.Parse(a.Substring(22, 4)) >= first_res_id && int.Parse(a.Substring(22, 4)) <= last_res_id).ToList();

            //this.residue_index = int.Parse(pdb_atom_line.Substring(22, 4));
            //this.i_code = pdb_atom_line[26];

            var output_pdb_file = $@"{pdb_out_folder}{pdb_id}{chain_id}_{first_res_id}_{last_res_id}.pdb";

            if (!File.Exists(output_pdb_file) || new FileInfo(output_pdb_file).Length == 0)
            {
                Directory.CreateDirectory(Path.GetDirectoryName(output_pdb_file));

                File.WriteAllLines(output_pdb_file, pdb_lines.ToList());
            }

            return output_pdb_file;
        }

        public static void run_dssp(string pdb_id, string pdb_folder = @"c:\betastrands_dataset\pdb\", string dssp_exe = @"c:\betastrands_dataset\pdb\dssp2.exe")
        {
            var pdb_file = $@"{pdb_folder}{pdb_id}.pdb";
            var dssp_file = $@"{pdb_folder}{pdb_id}.dssp";


            var start = new ProcessStartInfo
            {
                WorkingDirectory = Path.GetDirectoryName(dssp_exe) ?? "",
                FileName = dssp_exe,
                Arguments = $"-i \"{pdb_file}\" -o \"{dssp_file}\"",
                UseShellExecute = false,
                CreateNoWindow = false,
                RedirectStandardOutput = true,
                RedirectStandardError = true
            };

            using (var process = Process.Start(start))
            {
                if (process == null) return;

                Program.WriteLine(start.FileName + " " + start.Arguments);
                using (var reader = process.StandardOutput)
                {
                    var result = reader.ReadToEnd(); // Here is the result of StdOut(for example: print "test")
                    if (!string.IsNullOrWhiteSpace(result)) Program.WriteLine(result);

                    var stderr = process.StandardError.ReadToEnd(); // Here are the exceptions from our Python script
                    if (!string.IsNullOrWhiteSpace(stderr)) Program.WriteLine(stderr);
                    return;
                }
            }
        }

        public Atom(string pdb_id, string pdb_atom_line, int pdb_model_array_index = -1, int array_index = -1)
        {
            this.pdb_id = pdb_id;
            this.chain_id = pdb_atom_line[21];
            this.residue_index = int.Parse(pdb_atom_line.Substring(22, 4));
            this.i_code = pdb_atom_line[26];
            this.amino_acid = Atom.Aa3To1(pdb_atom_line.Substring(17, 3));
            this.array_index = array_index;
            this.model_index = pdb_model_array_index;

            this.serial_index = int.Parse(pdb_atom_line.Substring(7 - 1, (11 - 7) + 1));

            this.atom_type = pdb_atom_line.Substring(13 - 1, (16 - 13) + 1).Trim(); // is this CA / CB  etc ??
            this.element = pdb_atom_line.Length >= 77 - 1 ? pdb_atom_line.Substring(77 - 1, (78 - 77) + 1).Trim().FirstOrDefault() : default(char);
            this.X = double.Parse(pdb_atom_line.Substring(31 - 1, (38 - 31) + 1));
            this.Y = double.Parse(pdb_atom_line.Substring(39 - 1, (46 - 39) + 1));
            this.Z = double.Parse(pdb_atom_line.Substring(47 - 1, (54 - 47) + 1));
        }

        //public static bool use_ram_disk = true;
        public const string pdb_folder = /*$@"{(use_ram_disk?"r":"c")}*/ @"c:\betastrands_dataset\pdb\";


        public static List<(string pdb_id, int pdb_model_index, char chain_id, List<Atom> pdb_model_chain_atoms)> load_atoms_pdb(string pdb_id, bool select_first_icode = true, bool find_intramolecular = true, bool find_intermolecular = false, bool load_rsa_data = false, bool load_dssp_data = false, bool load_stride_data = false, string pdb_folder = Atom.pdb_folder)
        {
            var first_model_only = true;

            var load_ring_data = true;
            var load_mpsa_sec_struct_predictions = true;
            var load_blast_pssms = true;
            var load_iup_data = true;
            var load_sable = true;
            var load_ala_scan = true;
            var load_dna_binding_vars = true;

            pdb_id = Path.GetFileNameWithoutExtension(pdb_id);
            pdb_id = pdb_id.ToUpperInvariant();
            var pdb_id_simple = pdb_id.Substring(0, 4);
            var result = new List<(string pdb_id, int pdb_model_index, char chain_id, List<Atom> pdb_model_chain_atoms)>();
            var pdb_filename = $@"{pdb_folder}{pdb_id}.pdb";
            var pdb_lines = File.ReadAllLines(pdb_filename).ToList();
            var pdb_model_array_index = -1;


            while (pdb_lines.Count > 0)
            {
                // set model index
                pdb_model_array_index++;

                // find end of model index
                var endmdl_index = pdb_lines.FindIndex(a => a.StartsWith("ENDMDL"));

                // get all pdb lines for this model
                var pdb_model_lines = (endmdl_index > -1) ? pdb_lines.Take(endmdl_index + 1).ToList() : pdb_lines;

                // remove model from pdb_lines
                pdb_lines = pdb_lines.Skip(pdb_model_lines.Count).ToList();

                // filter model for ATOMS only
                pdb_model_lines = pdb_model_lines.Where(a => a.StartsWith("ATOM ")).ToList();

                // if no ATOMs, go to next model
                if (pdb_model_lines.Count == 0) continue;

                // convert ATOM records to ATOM classes
                var pdb_model_atoms = pdb_model_lines.Select((atom_struct, array_index) => new Atom(pdb_id, atom_struct, pdb_model_array_index, array_index)).ToList();

                // save the chain order, which may not be alphabetical, or such
                var chain_order_in_pdb_file = pdb_model_atoms.Select(atom => atom.chain_id).Distinct().Select((chain_id, array_index) => (array_index, chain_id)).ToList();
                pdb_model_atoms.ForEach(a => a.chain_order_in_pdb_file = chain_order_in_pdb_file);

                // re-order the atom classes by chain id, res_id, i_code
                pdb_model_atoms = pdb_model_atoms.OrderBy(a => a.chain_id).ThenBy(a => a.residue_index).ThenBy(a => a.i_code).ToList();

                // remove any additional i_codes
                if (select_first_icode)
                {
                    pdb_model_atoms = pdb_model_atoms.GroupBy(a => (a.chain_id, a.residue_index)).SelectMany(a => a.ToList().GroupBy(b => b.i_code).First().ToList()).ToList();
                }

                // make list of atoms for each amino acid
                foreach (var atom in pdb_model_atoms)
                {
                    if (atom.amino_acid_atoms != null) continue;
                    atom.amino_acid_atoms = pdb_model_atoms.Where(a => a.model_index == atom.model_index && string.Equals(a.pdb_id, atom.pdb_id, StringComparison.InvariantCultureIgnoreCase) && a.chain_id == atom.chain_id && a.residue_index == atom.residue_index && a.i_code == atom.i_code).ToList();
                    atom.amino_acid_atoms.ForEach(a => a.amino_acid_atoms = atom.amino_acid_atoms);
                }

                // get the main/master atoms
                var pdb_model_master_atoms = Atom.select_amino_acid_master_atoms(pdb_id, pdb_model_atoms);

                // set the master_atom for all atoms
                for (var master_atom_index = 0; master_atom_index < pdb_model_master_atoms.Count; master_atom_index++)
                {
                    var master_atom = pdb_model_master_atoms[master_atom_index];

                    master_atom.amino_acid_atoms.ForEach(a => a.amino_acid_master_atom = master_atom);
                }

                // load ring values
                if (load_ring_data) Atom.load_ring(pdb_id, pdb_model_atoms);

                // load dssp values
                if (load_dssp_data) Atom.load_dssp(pdb_id, pdb_model_atoms);

                // load stride values
                if (load_stride_data) Atom.load_stride(pdb_id, pdb_model_atoms);

                // load free-sasa values (rsa)
                if (load_rsa_data) Atom.load_rsa(pdb_id, pdb_model_atoms);


                // load psi-blast PSSMs

                if (load_blast_pssms) Atom.load_pssm(pdb_id, pdb_model_atoms);

                // load MPSA secondary structure predictions
                if (load_mpsa_sec_struct_predictions)
                {
                    Atom.load_mpsa_sec_struct_predictions(pdb_id, pdb_model_atoms);

                    //pdb_model_chain_atoms.ForEach(a => compare_dssp_to_mpsa(pdb_id, a.chain_id, a.pdb_model_chain_atoms));
                }


                // find intra/inter-molecular contacts/interactions
                if (find_intermolecular) Atom.find_intermolecular_contacts(pdb_id, pdb_model_atoms, (double) 5);
                if (find_intramolecular) Atom.find_intramolecular_contacts(pdb_id, pdb_model_atoms, (double) 5);


                // further split model atoms into chain atoms
                var pdb_model_chain_atoms = pdb_model_atoms.GroupBy(a => a.chain_id).Select(a => (chain_id: a.Key, pdb_model_chain_atoms: a.ToList())).ToList();
                var pdb_model_chain_master_atoms = pdb_model_chain_atoms.Select(a => (chain_id: a.chain_id, pdb_model_chain_master_atoms: select_amino_acid_master_atoms(pdb_id, a.pdb_model_chain_atoms))).ToList();

                foreach (var c in pdb_model_chain_master_atoms)
                {
                    // get the main/master atoms
                    //var pdb_model_chain_master_atoms = Atom.select_amino_acid_master_atoms(pdb_id, c.pdb_model_chain_atoms);


                    // set the master_atom for all atoms
                    for (var master_atom_index = 0; master_atom_index < c.pdb_model_chain_master_atoms.Count; master_atom_index++)
                    {
                        var master_atom = c.pdb_model_chain_master_atoms[master_atom_index];

                        master_atom.amino_acid_atoms.ForEach(a => a.amino_acid_master_atom = master_atom);
                        master_atom.amino_acid_atoms.ForEach(a => a.master_index = master_atom_index);
                    }
                }

                if (load_iup_data)
                {
                    foreach (var c in pdb_model_chain_master_atoms)
                    {
                        // get the main/ master atoms
                        //var pdb_model_chain_master_atoms = Atom.select_amino_acid_master_atoms(pdb_id, c.pdb_model_chain_atoms);

                        //var chain_seq = string.Join("", pdb_model_chain_master_atoms.Select(a => a.amino_acid).ToList());


                        var iup_data = iupred.load(pdb_id_simple, c.chain_id);

                        if (iup_data == null || iup_data.Count == 0)
                        {
                            continue;
                        }

                        for (var master_atom_index = 0; master_atom_index < c.pdb_model_chain_master_atoms.Count; master_atom_index++)
                        {
                            var master_atom = c.pdb_model_chain_master_atoms[master_atom_index];

                            master_atom.iup_entry = iup_data[master_atom_index];

                            foreach (var atom in c.pdb_model_chain_master_atoms[master_atom_index].amino_acid_atoms)
                            {
                                atom.iup_entry = master_atom.iup_entry;
                            }
                        }
                    }
                }

                if (load_dna_binding_vars)
                {
                    foreach (var c in pdb_model_chain_master_atoms)
                    {
                        load_dna_binding(pdb_id + c.chain_id, c.pdb_model_chain_master_atoms);
                    }
                }

                if (load_ala_scan)
                {
                    foreach (var c in pdb_model_chain_atoms)
                    {
                        // load foldx AlaScan values

                        // 1EV0A_Repair_AS.fxout
                        var alascan = foldx_caller.load_foldx_ala_scanning(pdb_id + c.chain_id + "_Repair", c.chain_id, null, false).data;

                        if (alascan != null && alascan.Count > 0)
                        {
                            foreach (var scan in alascan)
                            {
                                var atoms = c.pdb_model_chain_atoms.Where(a => string.Equals(a.pdb_id + a.chain_id + "_Repair", scan.pdb_id, StringComparison.InvariantCultureIgnoreCase) && a.chain_id == scan.chain_id && a.residue_index == scan.residue_index).ToList();
                                atoms.ForEach(a => a.foldx_monomer_ala_scan_ddg = scan.ddg);
                            }
                        }
                    }

                    // load foldx PositionScan values
                    // var load_pos_scan = true;
                    // if (load_pos_scan) Atom.load_foldx_position_scanning(pdb_model_atoms);
                }

                // load uniprot sequence (& still to do: also load mapping and save in atoms)
                //var chain_ids = pdb_model_chain_atoms.Select(a => a.chain_id).Distinct().ToList();
                foreach (var c in pdb_model_chain_atoms)
                    //    foreach (var chain_id in chain_ids)
                {
                    var uniprot_file = $@"c:\betastrands_dataset\uniprot\{pdb_id_simple}{c.chain_id}.fasta";

                    if (!File.Exists(uniprot_file) || new FileInfo(uniprot_file).Length == 0) continue;

                    var uniprot_sequence = string.Join("", File.ReadAllLines(uniprot_file).Where(a => !string.IsNullOrWhiteSpace(a) && !a.StartsWith(">")).ToList());

                    c.pdb_model_chain_atoms.ForEach(a => a.uniprot_sequence = uniprot_sequence);
                }

                // load sable
                if (load_sable)
                {
                    foreach (var c in pdb_model_chain_master_atoms)
                    {
                        var sable_file = $@"C:\betastrands_dataset\sable\{pdb_id_simple}{c.chain_id}.txt";
                        var sable_data = sable.load(sable_file);


                        // check sable sequences matches this pdb
                        for (var index = 0; index < c.pdb_model_chain_master_atoms.Count; index++)
                        {
                            c.pdb_model_chain_master_atoms[index].sable_item = sable_data[index];

                            if (c.pdb_model_chain_master_atoms[index].amino_acid != sable_data[index].amino_acid) throw new Exception();
                        }

                        // load sable data into atoms
                        for (var index = 0; index < c.pdb_model_chain_master_atoms.Count; index++)
                        {
                            var master_atom = c.pdb_model_chain_master_atoms[index];

                            master_atom.sable_item = sable_data[index];
                        }

                        // load ss as mpsa instance
                        var list = sable_data.Select(a => (a.seq_index, a.amino_acid, a.predicted_ss, a.prob_h, a.prob_e, a.prob_e)).ToList();
                        var mpsa_reader_sable = new mpsa_reader("sable", list);

                        for (var index = 0; index < c.pdb_model_chain_master_atoms.Count; index++)
                        {
                            var master_atom = c.pdb_model_chain_master_atoms[index];
                            if (master_atom.mpsa_entries == null) master_atom.mpsa_entries = new List<(string format, mpsa_reader.mpsa_line_entry)>();
                            var mpsa_matrix_line_entry = mpsa_reader_sable?.mpsa_matrix?.First /*OrDefault*/(a => a.index == index);


                            master_atom.mpsa_entries.Add((mpsa_reader_sable.format, mpsa_matrix_line_entry));

                            foreach (var atom in master_atom.amino_acid_atoms)
                            {
                                atom.mpsa_entries = master_atom.mpsa_entries;
                            }
                        }
                    }
                }

                // add split model/chain to returned result
                pdb_model_chain_atoms.ForEach(a => result.Add((pdb_id, pdb_model_array_index, a.chain_id, a.pdb_model_chain_atoms)));


                // set chain_atoms property
                pdb_model_chain_atoms.ForEach(a => a.pdb_model_chain_atoms.ForEach(b => b.chain_atoms = a.pdb_model_chain_atoms));

                // output data for comparison of DSSP & STRIDE data to MPSA secondary structure predictions
                //pdb_model_chain_atoms.ForEach(a => compare_dssp_to_mpsa(pdb_id, a.chain_id, a.pdb_model_chain_atoms));
            }

            if (first_model_only)
            {
                result = result.Where(a => a.pdb_model_index == 0).ToList();
            }

            foreach (var r in result)
            {
                var dssp_mpsa = get_dssp_and_mpsa_subsequences(r.pdb_id, r.chain_id, r.pdb_model_chain_atoms, (get_dssp_and_mpsa_subsequences_params)0b_1111_1111_1111);

                var ground_truths = dssp_mpsa;// dssp_mpsa.Where(a => a.format.StartsWith("monomer", StringComparison.InvariantCultureIgnoreCase) || a.format.StartsWith("multimer", StringComparison.InvariantCultureIgnoreCase)).ToList();

                var data = new List<string>();

                var header = "pdb,chain_id,format,prediction," + string.Join(",",ground_truths.Select(a=>a.format).ToList());

                data.Add(header);

                foreach (var a in dssp_mpsa)
                {
                    var x = $"{r.pdb_id},{r.chain_id},{a.format},{a.prediction}";

                    foreach (var g in ground_truths)
                    {
                        var score = 0;
                        for (var i = 0; i < g.prediction.Length; i++)
                        {
                            if (g.prediction[i] == a.prediction[i]) score++;
                        }

                        var p = (double) score / (double) g.prediction.Length;

                        x = x + "," + p;
                    }

                    data.Add(x);
                }

                data.Add("");

                lock (dssp_mpsa_lock)
                {
                    File.AppendAllLines($@"c:\betastrands_dataset\dssp_mpsa_protein_q3_data.csv", data);
                }
            }

            return result;
        }

        public static object dssp_mpsa_lock = new object();

        public static List<string> compare_dssp_to_mpsa_list = new List<string>();

        public static object compare_dssp_to_mpsa_lock = new object();

        public static void compare_dssp_to_mpsa(string pdb_id, char chain_id, List<Atom> param_atoms)
        {
            lock (compare_dssp_to_mpsa_lock)
            {
                if (compare_dssp_to_mpsa_list.Contains($"{pdb_id}{chain_id}")) return;
                compare_dssp_to_mpsa_list.Add($"{pdb_id}{chain_id}");

                var result1 = new List<(string filename, string pdb, string format, string data_category1, string data_category2, string[] data)>();
                var result2 = new List<(string filename, string pdb, string format, string data_category1, string data_category2, string[] data)>();

                var pdb = $"{pdb_id}{chain_id}";

                var master_atoms = Atom.select_amino_acid_master_atoms(null, param_atoms);

                var filename = "";
                var format = "pdb";
                //var blank = (filename, "", "", "", "", new string[] { });

                var aa_seq = string.Join("", master_atoms.Select(a => a.amino_acid).ToArray());
                var aa_seq_aa_types = aa_seq.Distinct().OrderBy(a => a).ToArray();

                var strand_interface_type_seq = string.Join("", master_atoms.Select(a => a.strand_interface_type).ToArray());
                var strand_interface_type_seq_aa_types = strand_interface_type_seq.Distinct().OrderBy(a => a).ToArray();

                var parallelism_seq = string.Join("", master_atoms.Select(a => a.parallelism).ToArray());
                var parallelism_seq_aa_types = parallelism_seq.Distinct().OrderBy(a => a).ToArray();

                var symmetry_mode_seq = string.Join("", master_atoms.Select(a => a.symmetry_mode).ToArray());
                var symmetry_mode_seq_aa_types = symmetry_mode_seq.Distinct().OrderBy(a => a).ToArray();

                var dssp_seq = string.Join("", master_atoms.Select(a => a.monomer_dssp).ToArray());
                var dssp_seq_ss_types = dssp_seq.Distinct().OrderBy(a => a).ToArray();

                var stride_seq = string.Join("", master_atoms.Select(a => a.monomer_stride).ToArray());
                var stride_seq_ss_types = stride_seq.Distinct().OrderBy(a => a).ToArray();

                var dssp3_seq = string.Join("", master_atoms.Select(a => a.monomer_dssp3).ToArray());
                var dssp3_seq_ss_types = dssp3_seq.Distinct().OrderBy(a => a).ToArray();

                var stride3_seq = string.Join("", master_atoms.Select(a => a.monomer_stride3).ToArray());
                var stride3_seq_ss_types = stride3_seq.Distinct().OrderBy(a => a).ToArray();


                // add sequences to result
                result1.Add((filename, pdb, format, "seq_index", "seq_index", aa_seq.Select((b, i) => (i + 1).ToString()).ToArray()));

                result1.Add((filename, pdb, format, "seq_aa", "seq_aa", aa_seq.Select(b => b.ToString()).ToArray()));
                result1.Add((filename, pdb, format, "seq_bsi", "seq_bsi", strand_interface_type_seq.Select(b => b.ToString()).ToArray()));
                result1.Add((filename, pdb, format, "seq_symmetry", "seq_symmetry", symmetry_mode_seq.Select(b => b.ToString()).ToArray()));
                result1.Add((filename, pdb, format, "seq_parallel", "seq_parallel", parallelism_seq.Select(b => b.ToString()).ToArray()));
                result1.Add((filename, pdb, format, "seq_dssp", "seq_dssp", dssp_seq.Select(b => b.ToString()).ToArray()));
                result1.Add((filename, pdb, format, "seq_dssp3", "seq_dssp3", dssp3_seq.Select(b => b.ToString()).ToArray()));
                result1.Add((filename, pdb, format, "seq_stride", "seq_stride", stride_seq.Select(b => b.ToString()).ToArray()));
                result1.Add((filename, pdb, format, "seq_stride3", "seq_stride3", stride3_seq.Select(b => b.ToString()).ToArray()));

                // add numeric values to result
                //aa_seq_aa_types.ToList().ForEach(aa => result.Add((filename, pdb, format, nameof(aa_seq), aa.ToString(), aa_seq.Select(a => aa == a ? 1.0 : 0.0).Select(b => b.ToString()).ToArray())));
                //strand_interface_type_seq_aa_types.Where(a => a != '_').ToList().ForEach(aa => result2.Add((filename, pdb, format, "seq_bsi", aa.ToString(), strand_interface_type_seq.Select(a => aa == a ? 1.0 : 0.0).Select(b => b.ToString()).ToArray())));
                //symmetry_mode_seq_aa_types.Where(a => a != '_').ToList().ForEach(aa => result2.Add((filename, pdb, format, "seq_symmetry", aa.ToString(), symmetry_mode_seq.Select(a => aa == a ? 1.0 : 0.0).Select(b => b.ToString()).ToArray())));                
                //parallelism_seq_aa_types.Where(a => a != '_').ToList().ForEach(aa => result2.Add((filename, pdb, format, "seq_parallel", aa.ToString(), parallelism_seq.Select(a => aa == a ? 1.0 : 0.0).Select(b => b.ToString()).ToArray())));                
                //dssp_seq_ss_types.ToList().ForEach(aa => result2.Add((filename, pdb, format, nameof(dssp_seq), aa.ToString(), dssp_seq.Select(a => aa == a ? 1.0 : 0.0).Select(b => b.ToString()).ToArray())));
                //dssp3_seq_ss_types.ToList().ForEach(aa => result2.Add((filename, pdb, format, nameof(dssp3_seq), aa.ToString(), dssp3_seq.Select(a => aa == a ? 1.0 : 0.0).Select(b => b.ToString()).ToArray())));
                //stride_seq_ss_types.ToList().ForEach(aa => result2.Add((filename, pdb, format, nameof(stride_seq), aa.ToString(), stride_seq.Select(a => aa == a ? 1.0 : 0.0).Select(b => b.ToString()).ToArray())));
                //stride3_seq_ss_types.ToList().ForEach(aa => result2.Add((filename, pdb, format, nameof(stride3_seq), aa.ToString(), stride3_seq.Select(a => aa == a ? 1.0 : 0.0).Select(b => b.ToString()).ToArray())));
                //var atoms_all_dssp_seq = string.Join("", atoms_all.Select(a => a.monomer_dssp).ToArray());
                //var atoms_all_dssp_seq_ss_types = atoms_all_dssp_seq.Distinct().OrderBy(a => a).ToArray();
                //var atoms_all_stride_seq = string.Join("", atoms_all.Select(a => a.monomer_stride).ToArray());
                //var atoms_all_stride_seq_ss_types = atoms_all_stride_seq.Distinct().OrderBy(a => a).ToArray();
                //var atoms_all_dssp3_seq = string.Join("", atoms_all.Select(a => a.monomer_dssp3).ToArray());
                //var atoms_all_dssp3_seq_ss_types = atoms_all_dssp3_seq.Distinct().OrderBy(a => a).ToArray();
                //var atoms_all_stride3_seq = string.Join("", atoms_all.Select(a => a.monomer_stride3).ToArray());
                //var atoms_all_stride3_seq_ss_types = atoms_all_stride3_seq.Distinct().OrderBy(a => a).ToArray();

                var mpsa_seqs = master_atoms.SelectMany(a => a.mpsa_entries).ToArray();

                var q3 = new List<(string pdb_id, char chain_id, string truth_format, string predictor_format, string aa_subset, char ss, double q3_value, double truth_total)>();

                q3.Add(("pdb_id",'c',"truth_format", "predictor_format", "aa_subset", 's', 0, 0));
                var av = new List<(string pdb_id, char chain_id, string predictor_format, string aa_subset, char ss, double average, double truth_total)>();


                // for each subset, what is the average predicted value
                av.Add(("pdb_id", 'c', "predictor_format", "aa_subset", 's', 0, 0));

                // problem(?): this doesn't separate the interfaces

                //const string  = ;

                foreach (var s in mpsa_seqs.GroupBy(a => a.format))
                {
                    format = s.Key;

                    {
                        // all
                        var atoms_all = master_atoms.Where(a => true).ToList();
                        var atoms_all_mpsa = atoms_all.SelectMany(a => a.mpsa_entries).Where(a => a.format == format).Select(a => a.mpsa_entry).ToArray();
                        var atoms_all_q3_dssp = "HEC*".Select(m => (ss: m, value: (double)atoms_all.Select((a, i) => ((m == '*' || atoms_all[i].monomer_dssp3 == m || atoms_all[i].monomer_dssp == m) && (atoms_all[i].monomer_dssp3 == atoms_all_mpsa[i].predicted_ss_code || atoms_all[i].monomer_dssp == atoms_all_mpsa[i].predicted_ss_code)) ? 1 : 0).Sum(), truth_total:(double)atoms_all.Count(b => m == '*' || b.monomer_dssp3 == m || b.monomer_dssp == m))).ToList();
                        atoms_all_q3_dssp.ForEach(a => q3.Add((pdb_id, chain_id, "dssp", format, nameof(atoms_all), a.ss, a.value, a.truth_total)));
                        var atoms_all_q3_stride = "HEC*".Select(m => (ss: m, value: (double)atoms_all.Select((a, i) => ((m == '*' || atoms_all[i].monomer_stride3 == m || atoms_all[i].monomer_stride == m) && (atoms_all[i].monomer_stride3 == atoms_all_mpsa[i].predicted_ss_code || atoms_all[i].monomer_stride == atoms_all_mpsa[i].predicted_ss_code)) ? 1 : 0).Sum(), truth_total:(double)atoms_all.Count(b => m == '*' || b.monomer_stride3 == m || b.monomer_stride == m))).ToList();
                        atoms_all_q3_stride.ForEach(a => q3.Add((pdb_id, chain_id, "stride", format, nameof(atoms_all), a.ss, a.value, a.truth_total)));
                        var atoms_all_mpsa_average = atoms_all_mpsa.SelectMany(a => a.line_prob_values).GroupBy(a => a.ss).Select(a => (ss: a.Key, av: a.Select(b => b.value).Average(), truth_total: a.Count())).ToList(); // each SS type, average probability
                        atoms_all_mpsa_average.ForEach(a => av.Add((pdb_id, chain_id, format, nameof(atoms_all), a.ss, a.av, a.truth_total)));
                        "HEC".ToList().ForEach(a=>av.Add((pdb_id, chain_id, "dssp", nameof(atoms_all),a, (double)atoms_all.Count(b=>b.monomer_dssp3==a || b.monomer_dssp==a) / (double)atoms_all.Count, atoms_all.Count))); //"HECT"
                    }

                    {
                        // non-interface
                        var atoms_non_interface = master_atoms.Where(a => a.strand_interface_type == '_').ToList();
                        var atoms_non_interface_mpsa = atoms_non_interface.SelectMany(a => a.mpsa_entries).Where(a => a.format == format).Select(a => a.mpsa_entry).ToArray();
                        var atoms_non_interface_q3_dssp = "HEC*".Select(m => (ss: m, value: (double)atoms_non_interface.Select((a, i) => ((m == '*' || atoms_non_interface[i].monomer_dssp3 == m || atoms_non_interface[i].monomer_dssp == m) && (atoms_non_interface[i].monomer_dssp3 == atoms_non_interface_mpsa[i].predicted_ss_code || atoms_non_interface[i].monomer_dssp == atoms_non_interface_mpsa[i].predicted_ss_code)) ? 1 : 0).Sum(), truth_total:(double)atoms_non_interface.Count(b => m == '*' || b.monomer_dssp3 == m || b.monomer_dssp == m))).ToList();
                        atoms_non_interface_q3_dssp.ForEach(a => q3.Add((pdb_id, chain_id, "dssp", format, nameof(atoms_non_interface), a.ss, a.value, a.truth_total)));
                        var atoms_non_interface_q3_stride = "HEC*".Select(m => (ss: m, value: (double)atoms_non_interface.Select((a, i) => ((m == '*' || atoms_non_interface[i].monomer_stride3 == m || atoms_non_interface[i].monomer_stride == m) && (atoms_non_interface[i].monomer_stride3 == atoms_non_interface_mpsa[i].predicted_ss_code || atoms_non_interface[i].monomer_stride == atoms_non_interface_mpsa[i].predicted_ss_code)) ? 1 : 0).Sum(), truth_total:(double)atoms_non_interface.Count(b => m == '*' || b.monomer_stride3 == m || b.monomer_stride == m))).ToList();
                        atoms_non_interface_q3_stride.ForEach(a => q3.Add((pdb_id, chain_id, "stride", format, nameof(atoms_non_interface), a.ss, a.value, a.truth_total)));
                        var atoms_non_interface_mpsa_average = atoms_non_interface_mpsa.SelectMany(a => a.line_prob_values).GroupBy(a => a.ss).Select(a => (ss: a.Key, av: a.Select(b => b.value).Average(), truth_total: a.Count())).ToList();
                        atoms_non_interface_mpsa_average.ForEach(a => av.Add((pdb_id, chain_id, format, nameof(atoms_non_interface), a.ss, a.av, a.truth_total)));
                    }

                    {
                        // dimrophic
                        var atoms_dimorphic_strand = master_atoms.Where((a, i) => a.strand_interface_type == 'S').ToList();
                        var atoms_dimorphic_strand_mpsa = atoms_dimorphic_strand.SelectMany(a => a.mpsa_entries).Where(a => a.format == format).Select(a => a.mpsa_entry).ToArray();
                        var atoms_dimorphic_strand_q3_dssp = "HEC*".Select(m => (ss: m, value: (double)atoms_dimorphic_strand.Select((a, i) => ((m == '*' || atoms_dimorphic_strand[i].monomer_dssp3 == m || atoms_dimorphic_strand[i].monomer_dssp == m) && (atoms_dimorphic_strand[i].monomer_dssp3 == atoms_dimorphic_strand_mpsa[i].predicted_ss_code || atoms_dimorphic_strand[i].monomer_dssp == atoms_dimorphic_strand_mpsa[i].predicted_ss_code)) ? 1 : 0).Sum(), truth_total:(double)atoms_dimorphic_strand.Count(b => m == '*' || b.monomer_dssp3 == m || b.monomer_dssp == m))).ToList();
                        atoms_dimorphic_strand_q3_dssp.ForEach(a => q3.Add((pdb_id, chain_id, "dssp", format, nameof(atoms_dimorphic_strand), a.ss, a.value, a.truth_total)));
                        var atoms_dimorphic_strand_q3_stride = "HEC*".Select(m => (ss: m, value: (double)atoms_dimorphic_strand.Select((a, i) => ((m == '*' || atoms_dimorphic_strand[i].monomer_stride3 == m || atoms_dimorphic_strand[i].monomer_stride == m) && (atoms_dimorphic_strand[i].monomer_stride3 == atoms_dimorphic_strand_mpsa[i].predicted_ss_code || atoms_dimorphic_strand[i].monomer_stride == atoms_dimorphic_strand_mpsa[i].predicted_ss_code)) ? 1 : 0).Sum(), truth_total:(double)atoms_dimorphic_strand.Count(b => m == '*' || b.monomer_stride3 == m || b.monomer_stride == m))).ToList();
                        atoms_dimorphic_strand_q3_stride.ForEach(a => q3.Add((pdb_id, chain_id, "stride", format, nameof(atoms_dimorphic_strand), a.ss, a.value, a.truth_total)));

                        var atoms_dimorphic_strand_mpsa_average = atoms_dimorphic_strand_mpsa.SelectMany(a => a.line_prob_values).GroupBy(a => a.ss).Select(a => (ss: a.Key, av: a.Select(b => b.value).Average(), truth_total: a.Count())).ToList();
                        atoms_dimorphic_strand_mpsa_average.ForEach(a => av.Add((pdb_id, chain_id, format, nameof(atoms_dimorphic_strand), a.ss, a.av, a.truth_total)));
                        var atoms_dimorphic_strand_flanking = master_atoms.Where((a, i) =>
                        {
                            // AA must not be strand interface (dimorphic/standard/hybrid)
                            if (a.strand_interface_type != '_') return false;

                            for (var j = 0; j <= 5; j++)
                            {
                                if (master_atoms.ElementAtOrDefault(i - j) != null && master_atoms[i - j].strand_interface_type == 'S') return true;
                                if (master_atoms.ElementAtOrDefault(i + j) != null && master_atoms[i + j].strand_interface_type == 'S') return true;
                            }

                            return false;
                        }).ToList();
                        var atoms_dimorphic_strand_flanking_mpsa = atoms_dimorphic_strand_flanking.SelectMany(a => a.mpsa_entries).Where(a => a.format == format).Select(a => a.mpsa_entry).ToArray();
                        var atoms_dimorphic_strand_flanking_q3_dssp = "HEC*".Select(m => (ss: m, value: (double)atoms_dimorphic_strand_flanking.Select((a, i) => ((m == '*' || atoms_dimorphic_strand_flanking[i].monomer_dssp3 == m || atoms_dimorphic_strand_flanking[i].monomer_dssp == m) && (atoms_dimorphic_strand_flanking[i].monomer_dssp3 == atoms_dimorphic_strand_flanking_mpsa[i].predicted_ss_code || atoms_dimorphic_strand_flanking[i].monomer_dssp == atoms_dimorphic_strand_flanking_mpsa[i].predicted_ss_code)) ? 1 : 0).Sum(), truth_total:(double)atoms_dimorphic_strand_flanking.Count(b => m == '*' || b.monomer_dssp3 == m || b.monomer_dssp == m))).ToList();
                        atoms_dimorphic_strand_flanking_q3_dssp.ForEach(a => q3.Add((pdb_id, chain_id, "dssp", format, nameof(atoms_dimorphic_strand_flanking), a.ss, a.value, a.truth_total)));
                        var atoms_dimorphic_strand_flanking_q3_stride = "HEC*".Select(m => (ss: m, value: (double)atoms_dimorphic_strand_flanking.Select((a, i) => ((m == '*' || atoms_dimorphic_strand_flanking[i].monomer_stride3 == m || atoms_dimorphic_strand_flanking[i].monomer_stride == m) && (atoms_dimorphic_strand_flanking[i].monomer_stride3 == atoms_dimorphic_strand_flanking_mpsa[i].predicted_ss_code || atoms_dimorphic_strand_flanking[i].monomer_stride == atoms_dimorphic_strand_flanking_mpsa[i].predicted_ss_code)) ? 1 : 0).Sum(), truth_total:(double)atoms_dimorphic_strand_flanking.Count(b => m == '*' || b.monomer_stride3 == m || b.monomer_stride == m))).ToList();
                        atoms_dimorphic_strand_flanking_q3_stride.ForEach(a => q3.Add((pdb_id, chain_id, "stride", format, nameof(atoms_dimorphic_strand_flanking), a.ss, a.value, a.truth_total)));
                        var atoms_dimorphic_strand_flanking_mpsa_average = atoms_dimorphic_strand_flanking_mpsa.SelectMany(a => a.line_prob_values).GroupBy(a => a.ss).Select(a => (ss: a.Key, av: a.Select(b => b.value).Average(), truth_total: a.Count())).ToList();
                        atoms_dimorphic_strand_flanking_mpsa_average.ForEach(a => av.Add((pdb_id, chain_id, format, nameof(atoms_dimorphic_strand_flanking), a.ss, a.av, a.truth_total)));
                    }

                    {
                        // standard
                        var atoms_standard_strand = master_atoms.Where((a, i) => a.strand_interface_type == 'M').ToList();
                        var atoms_standard_strand_mpsa = atoms_standard_strand.SelectMany(a => a.mpsa_entries).Where(a => a.format == format).Select(a => a.mpsa_entry).ToArray();
                        var atoms_standard_strand_q3_dssp = "HEC*".Select(m => (ss: m, value: (double) atoms_standard_strand.Select((a, i) => ((m == '*' || atoms_standard_strand[i].monomer_dssp3 == m || atoms_standard_strand[i].monomer_dssp == m) && (atoms_standard_strand[i].monomer_dssp3 == atoms_standard_strand_mpsa[i].predicted_ss_code || atoms_standard_strand[i].monomer_dssp == atoms_standard_strand_mpsa[i].predicted_ss_code)) ? 1 : 0).Sum(), truth_total: (double) atoms_standard_strand.Count(b => m == '*' || b.monomer_dssp3 == m || b.monomer_dssp == m))).ToList();
                        atoms_standard_strand_q3_dssp.ForEach(a => q3.Add((pdb_id, chain_id, "dssp", format, nameof(atoms_standard_strand), a.ss, a.value, a.truth_total)));
                        var atoms_standard_strand_q3_stride = "HEC*".Select(m => (ss: m, value: (double) atoms_standard_strand.Select((a, i) => ((m == '*' || atoms_standard_strand[i].monomer_stride3 == m || atoms_standard_strand[i].monomer_stride == m) && (atoms_standard_strand[i].monomer_stride3 == atoms_standard_strand_mpsa[i].predicted_ss_code || atoms_standard_strand[i].monomer_stride == atoms_standard_strand_mpsa[i].predicted_ss_code)) ? 1 : 0).Sum(), truth_total: (double) atoms_standard_strand.Count(b => m == '*' || b.monomer_stride3 == m || b.monomer_stride == m))).ToList();
                        atoms_standard_strand_q3_stride.ForEach(a => q3.Add((pdb_id, chain_id, "stride", format, nameof(atoms_standard_strand), a.ss, a.value, a.truth_total)));
                        var atoms_standard_strand_mpsa_average = atoms_standard_strand_mpsa.SelectMany(a => a.line_prob_values).GroupBy(a => a.ss).Select(a => (ss: a.Key, av: a.Select(b => b.value).Average(), truth_total: a.Count())).ToList();
                        atoms_standard_strand_mpsa_average.ForEach(a => av.Add((pdb_id, chain_id, format, nameof(atoms_standard_strand), a.ss, a.av, a.truth_total)));
                        var atoms_standard_strand_flanking = master_atoms.Where((a, i) =>
                        {
                            // AA must not be strand interface (dimorphic/standard/hybrid)
                            if (a.strand_interface_type != '_') return false;

                            for (var j = 0; j <= 5; j++)
                            {
                                if (master_atoms.ElementAtOrDefault(i - j) != null && master_atoms[i - j].strand_interface_type == 'M') return true;
                                if (master_atoms.ElementAtOrDefault(i + j) != null && master_atoms[i + j].strand_interface_type == 'M') return true;
                            }

                            return false;
                        }).ToList();
                        var atoms_standard_strand_flanking_mpsa = atoms_standard_strand_flanking.SelectMany(a => a.mpsa_entries).Where(a => a.format == format).Select(a => a.mpsa_entry).ToArray();
                        var atoms_standard_strand_flanking_q3_dssp = "HEC*".Select(m => (ss: m, value: (double) atoms_standard_strand_flanking.Select((a, i) => ((m == '*' || atoms_standard_strand_flanking[i].monomer_dssp3 == m || atoms_standard_strand_flanking[i].monomer_dssp == m) && (atoms_standard_strand_flanking[i].monomer_dssp3 == atoms_standard_strand_flanking_mpsa[i].predicted_ss_code || atoms_standard_strand_flanking[i].monomer_dssp == atoms_standard_strand_flanking_mpsa[i].predicted_ss_code)) ? 1 : 0).Sum(), truth_total: (double) atoms_standard_strand_flanking.Count(b => m == '*' || b.monomer_dssp3 == m || b.monomer_dssp == m))).ToList();
                        atoms_standard_strand_flanking_q3_dssp.ForEach(a => q3.Add((pdb_id, chain_id, "dssp", format, nameof(atoms_standard_strand_flanking), a.ss, a.value, a.truth_total)));
                        var atoms_standard_strand_flanking_q3_stride = "HEC*".Select(m => (ss: m, value: (double) atoms_standard_strand_flanking.Select((a, i) => ((m == '*' || atoms_standard_strand_flanking[i].monomer_stride3 == m || atoms_standard_strand_flanking[i].monomer_stride == m) && (atoms_standard_strand_flanking[i].monomer_stride3 == atoms_standard_strand_flanking_mpsa[i].predicted_ss_code || atoms_standard_strand_flanking[i].monomer_stride == atoms_standard_strand_flanking_mpsa[i].predicted_ss_code)) ? 1 : 0).Sum(), truth_total: (double) atoms_standard_strand_flanking.Count(b => m == '*' || b.monomer_stride3 == m || b.monomer_stride == m))).ToList();
                        atoms_standard_strand_flanking_q3_stride.ForEach(a => q3.Add((pdb_id, chain_id, "stride", format, nameof(atoms_standard_strand_flanking), a.ss, a.value, a.truth_total)));
                        var atoms_dimorphic_strand_flanking_mpsa_average = atoms_standard_strand_flanking_mpsa.SelectMany(a => a.line_prob_values).GroupBy(a => a.ss).Select(a => (ss: a.Key, av: a.Select(b => b.value).Average(), truth_total: a.Count())).ToList();
                        atoms_dimorphic_strand_flanking_mpsa_average.ForEach(a => av.Add((pdb_id, chain_id, format, nameof(atoms_standard_strand_flanking), a.ss, a.av, a.truth_total)));

                        //if (atoms_standard_strand.Count > 0)
                        //{
                        //    Console.WriteLine("");
                        //}
                    }


                    {
                        // hybrid
                        var atoms_hybrid_strand = master_atoms.Where((a, i) => a.strand_interface_type == 'H').ToList();
                        var atoms_hybrid_strand_mpsa = atoms_hybrid_strand.SelectMany(a => a.mpsa_entries).Where(a => a.format == format).Select(a => a.mpsa_entry).ToArray();
                        var atoms_hybrid_strand_q3_dssp = "HEC*".Select(m => (ss: m, value: (double)atoms_hybrid_strand.Select((a, i) => ((m == '*' || atoms_hybrid_strand[i].monomer_dssp3 == m || atoms_hybrid_strand[i].monomer_dssp == m) && (atoms_hybrid_strand[i].monomer_dssp3 == atoms_hybrid_strand_mpsa[i].predicted_ss_code || atoms_hybrid_strand[i].monomer_dssp == atoms_hybrid_strand_mpsa[i].predicted_ss_code)) ? 1 : 0).Sum(), truth_total:(double)atoms_hybrid_strand.Count(b => m == '*' || b.monomer_dssp3 == m || b.monomer_dssp == m))).ToList();
                        atoms_hybrid_strand_q3_dssp.ForEach(a => q3.Add((pdb_id, chain_id, "dssp", format, nameof(atoms_hybrid_strand), a.ss, a.value, a.truth_total)));
                        var atoms_hybrid_strand_q3_stride = "HEC*".Select(m => (ss: m, value: (double)atoms_hybrid_strand.Select((a, i) => ((m == '*' || atoms_hybrid_strand[i].monomer_stride3 == m || atoms_hybrid_strand[i].monomer_stride == m) && (atoms_hybrid_strand[i].monomer_stride3 == atoms_hybrid_strand_mpsa[i].predicted_ss_code || atoms_hybrid_strand[i].monomer_stride == atoms_hybrid_strand_mpsa[i].predicted_ss_code)) ? 1 : 0).Sum(), truth_total:(double)atoms_hybrid_strand.Count(b => m == '*' || b.monomer_stride3 == m || b.monomer_stride == m))).ToList();
                        atoms_hybrid_strand_q3_stride.ForEach(a => q3.Add((pdb_id, chain_id, "stride", format, nameof(atoms_hybrid_strand), a.ss, a.value, a.truth_total)));
                        var atoms_hybrid_strand_mpsa_average = atoms_hybrid_strand_mpsa.SelectMany(a => a.line_prob_values).GroupBy(a => a.ss).Select(a => (ss: a.Key, av: a.Select(b => b.value).Average(), truth_total: a.Count())).ToList();
                        atoms_hybrid_strand_mpsa_average.ForEach(a => av.Add((pdb_id, chain_id, format, nameof(atoms_hybrid_strand), a.ss, a.av, a.truth_total)));
                        var atoms_hybrid_strand_flanking = master_atoms.Where((a, i) =>
                        {
                            // AA must not be strand interface (dimorphic/standard/hybrid)
                            if (a.strand_interface_type != '_') return false;

                            for (var j = 0; j <= 5; j++)
                            {
                                if (master_atoms.ElementAtOrDefault(i - j) != null && master_atoms[i - j].strand_interface_type == 'H') return true;
                                if (master_atoms.ElementAtOrDefault(i + j) != null && master_atoms[i + j].strand_interface_type == 'H') return true;
                            }

                            return false;
                        }).ToList();
                        var atoms_hybrid_strand_flanking_mpsa = atoms_hybrid_strand_flanking.SelectMany(a => a.mpsa_entries).Where(a => a.format == format).Select(a => a.mpsa_entry).ToArray();
                        var atoms_hybrid_strand_flanking_q3_dssp = "HEC*".Select(m => (ss: m, value: (double)atoms_hybrid_strand_flanking.Select((a, i) => ((m == '*' || atoms_hybrid_strand_flanking[i].monomer_dssp3 == m || atoms_hybrid_strand_flanking[i].monomer_dssp == m) && (atoms_hybrid_strand_flanking[i].monomer_dssp3 == atoms_hybrid_strand_flanking_mpsa[i].predicted_ss_code || atoms_hybrid_strand_flanking[i].monomer_dssp == atoms_hybrid_strand_flanking_mpsa[i].predicted_ss_code)) ? 1 : 0).Sum(), truth_total:(double)atoms_hybrid_strand_flanking.Count(b => m == '*' || b.monomer_dssp3 == m || b.monomer_dssp == m))).ToList();
                        atoms_hybrid_strand_flanking_q3_dssp.ForEach(a => q3.Add((pdb_id, chain_id, "dssp", format, nameof(atoms_hybrid_strand_flanking), a.ss, a.value, a.truth_total)));
                        var atoms_hybrid_strand_flanking_q3_stride = "HEC*".Select(m => (ss: m, value: (double)atoms_hybrid_strand_flanking.Select((a, i) => ((m == '*' || atoms_hybrid_strand_flanking[i].monomer_stride3 == m || atoms_hybrid_strand_flanking[i].monomer_stride == m) && (atoms_hybrid_strand_flanking[i].monomer_stride3 == atoms_hybrid_strand_flanking_mpsa[i].predicted_ss_code || atoms_hybrid_strand_flanking[i].monomer_stride == atoms_hybrid_strand_flanking_mpsa[i].predicted_ss_code)) ? 1 : 0).Sum(), truth_total:(double)atoms_hybrid_strand_flanking.Count(b => m == '*' || b.monomer_stride3 == m || b.monomer_stride == m))).ToList();
                        atoms_hybrid_strand_flanking_q3_stride.ForEach(a => q3.Add((pdb_id, chain_id, "stride", format, nameof(atoms_hybrid_strand_flanking), a.ss, a.value, a.truth_total)));
                        var atoms_hybrid_strand_flanking_mpsa_average = atoms_hybrid_strand_flanking_mpsa.SelectMany(a => a.line_prob_values).GroupBy(a => a.ss).Select(a => (ss: a.Key, av: a.Select(b => b.value).Average(), truth_total: a.Count())).ToList();
                        atoms_hybrid_strand_flanking_mpsa_average.ForEach(a => av.Add((pdb_id, chain_id, format, nameof(atoms_hybrid_strand_flanking), a.ss, a.av, a.truth_total)));


                    }
                }


                //foreach (var s in mpsa_seqs.GroupBy(a => a.format))
                //{
                //    result1.Add(blank);

                //    format = s.Key;
                //    var x = s.ToList();

                //    // get the predicted SS sequence
                //    var predictor_seq = string.Join("", x.Select(a => a.mpsa_entry.predicted_ss_code).ToArray());
                //    var predictor_seq_ss_types = predictor_seq.Union(x.First().mpsa_entry.ss_column_headers).Distinct().OrderBy(a => a).ToArray();
                //    result1.Add((filename, pdb, format, $"seq_{format}", "seq_predicted", predictor_seq.Select(b => b.ToString()).ToArray()));

                //    //predictor_seq_ss_types.ToList().ForEach(aa => result2.Add((filename, pdb, format, "ss_binary", aa.ToString(), predictor_seq.Select(a => aa == a ? 1.0 : 0.0).Select(b => b.ToString()).ToArray())));

                //    predictor_seq_ss_types.ToList().ForEach(aa => result2.Add((filename, pdb, format, "ss_probability", aa.ToString(), x.Select(a => a.mpsa_entry.line_prob_values.FirstOrDefault(b => b.ss == aa).value).Select(b => b.ToString()).ToArray())));

                //    // distance between SS1 and SS1 for same amino acid on the same predictor - abs

                //    var distances = new List<(int index, double value)>();


                //    predictor_seq_ss_types.ToList().ForEach(aa1 => predictor_seq_ss_types.ToList().ForEach(aa2 =>
                //    {
                //        if (aa1 <= aa2) return;


                //        var d = x.Select((a, i) =>
                //        {
                //            var e = Math.Abs(a.mpsa_entry.line_prob_values.FirstOrDefault(b => b.ss == aa1).value - a.mpsa_entry.line_prob_values.FirstOrDefault(b => b.ss == aa2).value);

                //            distances.Add((i, e));

                //            return e;
                //        }).Select(b => b.ToString()).ToArray();


                //        result2.Add((filename, pdb, format, $"ss_probability_distance_abs", $"{aa1}_{aa2}", d));
                //    }));

                //    // 1. dimorphics interface average e%, h%, c%
                //    // 2. strand interface average e%, h%, c%
                //    // 3. protein average
                //    // 4. protein excl. strand interfaces average


                //    result2.Add((filename, pdb, format, $"ss_probability_distance_abs", $"average_all", distances.GroupBy(a => a.index).Select(a => a.Select(b => b.value).Average()).Select(b => b.ToString()).ToArray()));


                // non-abs
                //predictor_seq_ss_types.ToList().ForEach(aa1 => predictor_seq_ss_types.ToList().ForEach(aa2 =>
                //{
                //    //if (aa1 <= aa2) return;
                //
                //    result2.Add((filename, pdb, format, $"ss_probability_distance", $"{aa1}_{aa2}",
                //        x.Select(a =>
                //                (a.mpsa_entry.line_prob_values.FirstOrDefault(b => b.ss == aa1).value -
                //                 a.mpsa_entry.line_prob_values.FirstOrDefault(b => b.ss == aa2).value))
                //            .Select(b => b.ToString()).ToArray()));
                //}));


                //var score_dssp = 0d;
                //var score_stride = 0d;
                //for (var i = 0; i < dssp_seq.Length; i++)
                //{
                //    if (dssp3_seq[i] == predictor_seq[i] || dssp_seq[i] == predictor_seq[i]) score_dssp++;
                //    if (stride3_seq[i] == predictor_seq[i] || stride_seq[i] == predictor_seq[i]) score_stride++;

                //}

                //var dssp_pct = (double)score_dssp / (double)dssp_seq.Length;
                //var stride_pct = (double)score_stride / (double)dssp_seq.Length;

                //result1.Add((filename, pdb, format, "q3_dssp", "q3_dssp", new string[] { dssp_pct.ToString() }));
                //result1.Add((filename, pdb, format, "q3_stride", "q3_stride", new string[] { stride_pct.ToString() }));

                //    result2.Add(blank);
                //}

                //result1.Add(blank);
                //result2.Add(blank);

                //result1.GroupBy(a => a.filename).ToList().ForEach(a => File.AppendAllLines($@"C:\betastrands_dataset\dssp_vs_mpsa\ss_predictor_comparison_{a.Key}.csv", a.Select(b => $"{b.pdb},{b.format},{b.data_category1},{b.data_category2},{string.Join(",", b.data)}").ToList()));
                //result2.GroupBy(a => a.filename).ToList().ForEach(a => File.AppendAllLines($@"C:\betastrands_dataset\dssp_vs_mpsa\ss_predictor_comparison_{a.Key}.csv", a.Select(b => $"{b.pdb},{b.format},{b.data_category1},{b.data_category2},{string.Join(",", b.data)}").ToList()));

                var f_q3 = $@"C:\betastrands_dataset\dssp_vs_mpsa\ss_predictor_comparison_{nameof(q3)}.csv";
                var r_q3 = q3.Select(a => $@"{a.pdb_id},{a.chain_id},{a.truth_format},{a.predictor_format},{a.aa_subset},{a.ss},{a.q3_value},{a.truth_total}").ToList();
                File.AppendAllLines(f_q3, r_q3);

                var f_av = $@"C:\betastrands_dataset\dssp_vs_mpsa\ss_predictor_comparison_{nameof(av)}.csv";
                var r_av = av.Select(a => $@"{a.pdb_id},{a.chain_id},{a.predictor_format},{a.aa_subset},{a.ss},{a.average},{a.truth_total}").ToList();
                File.AppendAllLines(f_av, r_av);
            }
        }

        [Flags]
        public enum get_dssp_and_mpsa_subsequences_params
        {
            none =               0b_0000_0000_0000,

            aa_seq =             0b_0000_0000_0001,

            monomer_dssp_seq =   0b_0000_0000_0010,
            monomer_stride_seq = 0b_0000_0000_0100,
            monomer_dssp3_seq =  0b_0000_0000_1000,
            monomer_stride3_seq= 0b_0000_0001_0000,
                                 
            multimer_dssp_seq=   0b_0000_0010_0000,
            multimer_stride_seq= 0b_0000_0100_0000,
            multimer_dssp3_seq=  0b_0000_1000_0000,
            multimer_stride3_seq=0b_0001_0000_0000,
        }

        public static List<(string format, string prediction)> get_dssp_and_mpsa_subsequences(string pdb_id, char chain_id, List<Atom> atoms, get_dssp_and_mpsa_subsequences_params get_dssp_and_mpsa_subsequences_params = get_dssp_and_mpsa_subsequences_params.none)
        {
            var result = new List<(string format, string prediction)>();

            //result.Add("DSSP vs MPSA: " + pdb_id + " " + chain_id);

            var master_atoms = Atom.select_amino_acid_master_atoms(null, atoms);

            var aa_seq = string.Join("", master_atoms.Select(a => a.amino_acid).ToList());

            var monomer_dssp_seq = string.Join("", master_atoms.Select(a => a.monomer_dssp).ToList());
            var monomer_stride_seq = string.Join("", master_atoms.Select(a => a.monomer_stride).ToList());
            var monomer_dssp3_seq = string.Join("", master_atoms.Select(a => a.monomer_dssp3).ToList());
            var monomer_stride3_seq = string.Join("", master_atoms.Select(a => a.monomer_stride3).ToList());

            var multimer_dssp_seq = string.Join("", master_atoms.Select(a => a.multimer_dssp).ToList());
            var multimer_stride_seq = string.Join("", master_atoms.Select(a => a.multimer_stride).ToList());
            var multimer_dssp3_seq = string.Join("", master_atoms.Select(a => a.multimer_dssp3).ToList());
            var multimer_stride3_seq = string.Join("", master_atoms.Select(a => a.multimer_stride3).ToList());

            if (get_dssp_and_mpsa_subsequences_params.HasFlag(get_dssp_and_mpsa_subsequences_params.aa_seq)) { result.Add((nameof(aa_seq), aa_seq)); }

            if (get_dssp_and_mpsa_subsequences_params.HasFlag(get_dssp_and_mpsa_subsequences_params.monomer_dssp_seq)) { result.Add((nameof(monomer_dssp_seq), monomer_dssp_seq)); }
            if (get_dssp_and_mpsa_subsequences_params.HasFlag(get_dssp_and_mpsa_subsequences_params.monomer_stride_seq)) { result.Add((nameof(monomer_stride_seq), monomer_stride_seq)); }
            if (get_dssp_and_mpsa_subsequences_params.HasFlag(get_dssp_and_mpsa_subsequences_params.monomer_dssp3_seq)) { result.Add((nameof(monomer_dssp3_seq), monomer_dssp3_seq)); }
            if (get_dssp_and_mpsa_subsequences_params.HasFlag(get_dssp_and_mpsa_subsequences_params.monomer_stride3_seq)) { result.Add((nameof(monomer_stride3_seq), monomer_stride3_seq)); }

            if (get_dssp_and_mpsa_subsequences_params.HasFlag(get_dssp_and_mpsa_subsequences_params.multimer_dssp_seq)) { result.Add((nameof(multimer_dssp_seq), multimer_dssp_seq)); }
            if (get_dssp_and_mpsa_subsequences_params.HasFlag(get_dssp_and_mpsa_subsequences_params.multimer_stride_seq)) { result.Add((nameof(multimer_stride_seq), multimer_stride_seq)); }
            if (get_dssp_and_mpsa_subsequences_params.HasFlag(get_dssp_and_mpsa_subsequences_params.multimer_dssp3_seq)) { result.Add((nameof(multimer_dssp3_seq), multimer_dssp3_seq)); }
            if (get_dssp_and_mpsa_subsequences_params.HasFlag(get_dssp_and_mpsa_subsequences_params.multimer_stride3_seq)) { result.Add((nameof(multimer_stride3_seq), multimer_stride3_seq)); }


            var mpsa_seqs = master_atoms.SelectMany(a => a.mpsa_entries).ToList();

            foreach (var s in mpsa_seqs.GroupBy(a => a.format))
            {
                var format = s.Key;
                var x = s.ToList();
                var mpsa_seq = string.Join("", x.Select(a => a.mpsa_entry.predicted_ss_code).ToList());

                result.Add((format, mpsa_seq));
            }

            var con_ss_seq = string.Join("", master_atoms.Select(a =>
            {
                var y = a.mpsa_entries.SelectMany(b => b.mpsa_entry.line_prob_values).ToList();

                var prob_h = y.Where(b => b.ss == 'H').Select(b => b.value).DefaultIfEmpty(0).Average();
                var prob_e = y.Where(b => b.ss == 'E').Select(b => b.value).DefaultIfEmpty(0).Average();
                var prob_c = y.Where(b => b.ss == 'C').Select(b => b.value).DefaultIfEmpty(0).Average();
                var prob_t = y.Where(b => b.ss == 'T').Select(b => b.value).DefaultIfEmpty(0).Average();
                //var prob_f = y.Where(b => b.ss == 'F').Select(b => b.value).DefaultIfEmpty(0).Average();

                var values = new List<(char ss, double value)>();
                values.Add(('C', prob_c));
                values.Add(('E', prob_e));
                values.Add(('H', prob_h));
                values.Add(('T', prob_t));
                //values.Add(('F', prob_f));
                values = values.OrderByDescending(c => c.value).ToList();

                var all_zero = values.All(d => d.value == 0);

                return all_zero ? 'C' : values.First().ss;
            }).ToList());

            result.Add(("consensus", con_ss_seq));

            return result;
        }


        public static List<Atom> select_amino_acid_master_atoms(string pdb_id, List<Atom> atoms)
        {
            atoms = atoms.Where(a => pdb_id == null || string.Equals(pdb_id, a.pdb_id, StringComparison.InvariantCultureIgnoreCase)).Distinct().ToList();

            atoms = atoms.OrderBy(a => a.chain_id).ThenBy(a => a.residue_index).ThenBy(a => a.i_code).ToList();

            atoms = atoms.GroupBy(atom => (atom.pdb_id, atom.chain_id, atom.residue_index /*, atom.i_code*/)).Select(grouped_atoms =>
            {
                var ca = grouped_atoms.FirstOrDefault(b => b.atom_type == "CA");
                if (ca != null) return ca;

                var n = grouped_atoms.FirstOrDefault(b => b.atom_type == "N");
                if (n != null) return n;

                var c = grouped_atoms.FirstOrDefault(b => b.atom_type == "C");
                if (c != null) return c;

                var o = grouped_atoms.FirstOrDefault(b => b.atom_type == "O");
                if (o != null) return o;

                var cb = grouped_atoms.FirstOrDefault(b => b.atom_type == "CB");
                if (cb != null) return cb;

                return grouped_atoms.FirstOrDefault(atom => atom != null);
            }).ToList();

            return atoms;
        }

        public static List<Atom> select_amino_acid_master_atoms_ignore_chain(string pdb_id, List<Atom> atoms)
        {
            atoms = atoms.Where(a => pdb_id == null || string.Equals(pdb_id, a.pdb_id, StringComparison.InvariantCultureIgnoreCase)).Distinct().ToList();

            atoms = atoms.OrderBy(a => a.residue_index).ThenBy(a => a.i_code).ToList();

            atoms = atoms.GroupBy(atom => (atom.pdb_id, atom.residue_index /*, atom.i_code*/)).Select(grouped_atoms =>
            {
                var ca = grouped_atoms.FirstOrDefault(b => b.atom_type == "CA");
                if (ca != null) return ca;

                var n = grouped_atoms.FirstOrDefault(b => b.atom_type == "N");
                if (n != null) return n;

                var c = grouped_atoms.FirstOrDefault(b => b.atom_type == "C");
                if (c != null) return c;

                var o = grouped_atoms.FirstOrDefault(b => b.atom_type == "O");
                if (o != null) return o;

                var cb = grouped_atoms.FirstOrDefault(b => b.atom_type == "CB");
                if (cb != null) return cb;

                return grouped_atoms.FirstOrDefault(atom => atom != null);
            }).ToList();

            return atoms;
        }

        public static List<(string pdb_id, char chain_id, List<Atom> atom_list)> group_atoms_by_pdb_id_and_chain_id(string pdb_id, List<Atom> atoms)
        {
            atoms = atoms.Distinct().ToList();
            atoms = atoms.OrderBy(a => a.chain_id).ThenBy(a => a.residue_index).ThenBy(a => a.i_code) /*.ThenBy(a => a.SerialIndex)*/.ToList();

            var seq = atoms.Where(a => string.Equals(a.pdb_id, pdb_id, StringComparison.InvariantCultureIgnoreCase)).GroupBy(a => (pdb_id: a.pdb_id.ToUpperInvariant(), chain_id: a.chain_id)).Select(a => (pdb_id: a.Key.pdb_id, chain_id: a.Key.chain_id, atom_list: a.ToList())).ToList();
            return seq;
        }

        public static List<(string pdb_id, char chain_id, List<Atom> atom_list)> group_master_atoms_by_pdb_id_and_chain_id(string pdb_id, List<Atom> atoms)
        {
            pdb_id = Path.GetFileNameWithoutExtension(pdb_id);


            atoms = select_amino_acid_master_atoms(pdb_id, atoms);
            var seq = atoms.GroupBy(a => (pdb_id: a.pdb_id.ToUpperInvariant(), chain_id: a.chain_id)).Select(a => (pdb_id: a.Key.pdb_id, chain_id: a.Key.chain_id, atom_list: a.ToList())).ToList();
            return seq;
        }

        public static List<(string pdb_id, char chain_id, string aa_sequence)> amino_acid_sequence(string pdb_id, List<Atom> atoms)
        {
            pdb_id = Path.GetFileNameWithoutExtension(pdb_id);

            var seq = group_master_atoms_by_pdb_id_and_chain_id(pdb_id, atoms).Select(a => (pdb_id: a.pdb_id, chain_id: a.chain_id, aa_sequence: String.Join("", a.atom_list.Select(b => b.amino_acid).ToList()))).ToList();
            return seq;
        }

        public enum structure_oligomisation
        {
            monomer,
            multimer
        }

        public static List<(string pdb_id, char chain_id, string ss_sequence)> ss_sequence(string pdb_id, List<Atom> atoms, structure_oligomisation structure_oligomisation, feature_calcs.ss_types ss_type = feature_calcs.ss_types.DSSP)
        {
            pdb_id = Path.GetFileNameWithoutExtension(pdb_id);


            return group_master_atoms_by_pdb_id_and_chain_id(pdb_id, atoms).Select(a => (pdb_id: a.pdb_id, chain_id: a.chain_id, ss_sequence: String.Join("", a.atom_list.Select(b =>
            {
                if (ss_type == feature_calcs.ss_types.DSSP) return structure_oligomisation == structure_oligomisation.multimer ? b.multimer_dssp : b.monomer_dssp;
                if (ss_type == feature_calcs.ss_types.DSSP3) return structure_oligomisation == structure_oligomisation.multimer ? b.multimer_dssp3 : b.monomer_dssp3;

                if (ss_type == feature_calcs.ss_types.STRIDE) return structure_oligomisation == structure_oligomisation.multimer ? b.multimer_stride : b.monomer_stride;
                if (ss_type == feature_calcs.ss_types.STRIDE3) return structure_oligomisation == structure_oligomisation.multimer ? b.multimer_stride3 : b.monomer_stride3;
                return ' ';
            }).ToList()))).ToList();
        }

        private static void load_mpsa_sec_struct_predictions(string pdb_id, List<Atom> pdb_model_atoms)
        {
            //this.mpsa_entries = new mpsa_reader("");

            var pdb_id_simple = Path.GetFileNameWithoutExtension(pdb_id).Substring(0, 4);

            var master_atoms = Atom.select_amino_acid_master_atoms(pdb_id, pdb_model_atoms);


            var lookup_sequence_table = File.ReadAllLines(@"C:\betastrands_dataset\betastrands_dataset_sequences.txt").ToList();


            var sequences = Atom.amino_acid_sequence(pdb_id, master_atoms);

            foreach (var seq in sequences)
            {
                var chain_id = seq.chain_id;

                var chain_master_atoms = master_atoms.Where(a => a.chain_id == chain_id).ToList();
                var sequence = seq.aa_sequence;

                //var line_index_match = lookup_sequence_table.IndexOf(sequence);

                //if (line_index_match < 0) throw new Exception();

                if (pdb_id_simple == "1EGP")
                {
                }

                var format_list = mpsa_reader.secondary_structure_codes.Select(a => a.format).ToList();

                //var mpsa_files_wc = $@"C:\betastrands_dataset\ss_meta\{line_index_match}.*";
                //var mpsa_files = Directory.GetFiles(Path.GetDirectoryName(mpsa_files_wc), Path.GetFileName(mpsa_files_wc));

                var mpsa_files = format_list.Select(a => $@"C:\betastrands_dataset\ss_meta\{pdb_id_simple}{chain_id}.{a}").ToList();

                foreach (var mpsa_file in mpsa_files)
                {
                    var reader = new mpsa_reader(mpsa_file, chain_master_atoms.Select(a => a.amino_acid).ToList());


                    for (var index = 0; index < chain_master_atoms.Count; index++)
                    {
                        var master_atom = chain_master_atoms[index];

                        if (master_atom.mpsa_entries == null) master_atom.mpsa_entries = new List<(string format, mpsa_reader.mpsa_line_entry)>();


                        var mpsa_matrix_line_entry = reader?.mpsa_matrix?.First /*OrDefault*/(a => a.index == index);

                        //if (mpsa_matrix == null || mpsa_matrix.line_prob_values==null || mpsa_matrix.line_prob_values.Count==0)
                        //{
                        //    mpsa_matrix = new mpsa_reader.mpsa_line_entry()
                        //    {   reader = reader,
                        //        index = index,
                        //        predicted_ss_code = ' ',
                        //        amino_acid = master_atom.amino_acid,
                        //        line_prob_values = reader.ss_column_headers.Select(a=>(a, master_atom.amino_acid, 0d)).ToList(),
                        //        ss_column_headers = reader.ss_column_headers
                        //    };
                        //}

                        master_atom.mpsa_entries.Add((reader.format, mpsa_matrix_line_entry));

                        foreach (var atom in master_atom.amino_acid_atoms)
                        {
                            atom.mpsa_entries = master_atom.mpsa_entries;
                        }
                    }

                    //}
                }

                // build consensus
            }
        }


        public static List<string> pssm_database_names = new List<string>();

        public static void load_dna_binding(string pdb_id, List<Atom> pdb_model_atoms)
        {
            var files = Directory.GetFiles(@"C:\betastrands_dataset\stackdppred_results\", $"{pdb_id}.txt", SearchOption.AllDirectories);

            foreach (var file in files)
            {
                var data = File.ReadAllLines(file).First();
                var data2 = string.Join("",data.Where(a => !"{}'':,".Contains(a)).ToList()).Split().ToList();

                var non_binding_prob = double.Parse(data2[data2.IndexOf("non-binding_prob") + 1]);
                var binding_prob = double.Parse(data2[data2.IndexOf("binding_prob") + 1]);

                var db = Path.GetDirectoryName(file).Split(new char[] { '\\', '/' },StringSplitOptions.RemoveEmptyEntries).Last().Split('_').Last();

                foreach (var atom in pdb_model_atoms)
                {
                    if (string.Equals(db, "nr", StringComparison.InvariantCultureIgnoreCase)) atom.chain_dna_binding_prob_nr = binding_prob;
                    if (string.Equals(db, "swissprot", StringComparison.InvariantCultureIgnoreCase)) atom.chain_dna_binding_prob_swissprot = binding_prob;
                    if (string.Equals(db, "uniref90", StringComparison.InvariantCultureIgnoreCase)) atom.chain_dna_binding_prob_uniref90 = binding_prob;
                }

                //Console.WriteLine();
                
            }
        }

        public static void load_pssm(string pdb_id, List<Atom> pdb_model_atoms)
        {
            var pdb_id_simple = Path.GetFileNameWithoutExtension(pdb_id).Substring(0, 4);

            var master_atoms = Atom.select_amino_acid_master_atoms(pdb_id, pdb_model_atoms);

            //var lookup_sequence_table = File.ReadAllLines(@"C:\betastrands_dataset\betastrands_dataset_sequences.txt").ToList();

            var sequences = Atom.amino_acid_sequence(pdb_id, master_atoms);

            foreach (var seq in sequences)
            {
                var chain_id = seq.chain_id;

                var chain_master_atoms = master_atoms.Where(a => a.chain_id == chain_id).ToList();
                var sequence = seq.aa_sequence;

                //var line_index_match = lookup_sequence_table.IndexOf(sequence);

                //if (line_index_match < 0) throw new Exception();

                var pssm_folders = Directory.GetDirectories($@"C:\betastrands_dataset\", "blast_pssm_*");
                var pssm_database_names = pssm_folders.Select(a => a.Split(new char[] {'\\', '/'}, StringSplitOptions.RemoveEmptyEntries).Last()).ToList();


                //var pssm_files = new string[]
                //    {
                //        $@"C:\betastrands_dataset\blast_pssm_swissprot_local\{line_index_match}.pssm",
                //        $@"C:\betastrands_dataset\blast_pssm_nr_remote\{line_index_match}.pssm"
                //    };

                var pssm_files = pssm_folders.Select(a => Path.Combine(a, $"{pdb_id_simple}{chain_id}.pssm")).ToArray();

                foreach (var pssm_file in pssm_files)
                {
                    var pssm_database_name = Path.GetDirectoryName(pssm_file).Split(new char[] {'\\', '/'}, StringSplitOptions.RemoveEmptyEntries).Last();

                    //if (File.Exists(pssm_file) && new FileInfo(pssm_file).Length > 0)
                    {
                        var pssm_matrix_unnormalised = pssm.load_psi_blast_pssm(pssm_file);
                        var pssm_matrix_normalised = pssm.normalise_pssm(pssm_matrix_unnormalised);

                        for (var index = 0; index < chain_master_atoms.Count; index++)
                        {
                            //var amino_acid_pssm = pssm_match.First(a => a.matrix_row_index == index);
                            var amino_acid_pssm_unnormalised = pssm_matrix_unnormalised.Where(a => a.matrix_row_index == index).ToList();
                            var amino_acid_pssm_normalised = pssm_matrix_normalised.Where(a => a.matrix_row_index == index).ToList();

                            if (amino_acid_pssm_unnormalised.Count != 20)
                            {
                                //throw new Exception();
                                amino_acid_pssm_unnormalised = new List<pssm.pssm_entry>();
                                amino_acid_pssm_normalised = new List<pssm.pssm_entry>();
                            }

                            var master_atom = chain_master_atoms[index];

                            if (master_atom.amino_acid_pssm_unnormalised == null) master_atom.amino_acid_pssm_unnormalised = new List<(string database, List<pssm.pssm_entry> pssm_entries)>();
                            if (master_atom.amino_acid_pssm_normalised == null) master_atom.amino_acid_pssm_normalised = new List<(string database, List<pssm.pssm_entry> pssm_entries)>();

                            master_atom.amino_acid_pssm_unnormalised.Add((pssm_database_name, amino_acid_pssm_unnormalised));
                            master_atom.amino_acid_pssm_normalised.Add((pssm_database_name, amino_acid_pssm_normalised));

                            foreach (var a in master_atom.amino_acid_atoms)
                            {
                                a.amino_acid_pssm_unnormalised = master_atom.amino_acid_pssm_unnormalised;
                                a.amino_acid_pssm_normalised = master_atom.amino_acid_pssm_normalised;
                            }
                        }
                    }
                }
            }
        }

        public static void load_ring(string pdb_id, List<Atom> atoms)
        {
            //var edge_files = Directory.GetFiles(@"C:\Users\aaron\Desktop\ring_pdb\pdb\", "*.edges");
            //var node_files = Directory.GetFiles(@"C:\Users\aaron\Desktop\ring_pdb\pdb\", "*.nodes");

            pdb_id = Path.GetFileNameWithoutExtension(pdb_id);

            var ring_folder = @"C:\betastrands_dataset\ring_pdb\pdb\";

            // Note: careful from the added '_Repair' in the filename

            var chains = atoms.Select(a => a.chain_id).Distinct().ToList();

            foreach (var chain in chains)
            {
                var ring_monomer_edge_filename = Path.Combine(ring_folder, pdb_id + chain + "_Repair.pdb.edges");
                var ring_monomer_node_filename = Path.Combine(ring_folder, pdb_id + chain + "_Repair.pdb.nodes");

                if (!File.Exists(ring_monomer_edge_filename) || new FileInfo(ring_monomer_edge_filename).Length == 0) Console.WriteLine("Warning: Ring Edge File is missing or empty: " + ring_monomer_edge_filename);
                if (!File.Exists(ring_monomer_node_filename) || new FileInfo(ring_monomer_node_filename).Length == 0) Console.WriteLine("Warning: Ring Node File is missing or empty: " + ring_monomer_edge_filename);

                var ring_monomer_edges = ring.ring_edge.load(ring_monomer_edge_filename);

                var ring_monomer_nodes = ring.ring_node.load(ring_monomer_node_filename);

                foreach (var ring_monomer_edge in ring_monomer_edges)
                {
                    var edge_atoms = atoms.Where(a => (ring_monomer_edge.NodeId1.chain == a.chain_id && ring_monomer_edge.NodeId1.res_id == a.residue_index && ring_monomer_edge.NodeId1.icode == a.i_code && ring_monomer_edge.NodeId1.amino_acid1 == a.amino_acid)
                        //||
                        //(ring_monomer_edge.NodeId2.chain == a.chain_id &&
                        // ring_monomer_edge.NodeId2.res_id == a.residue_index &&
                        // ring_monomer_edge.NodeId2.icode == a.i_code &&
                        // ring_monomer_edge.NodeId2.amino_acid1 == a.amino_acid)
                    ).ToList();

                    foreach (var a in edge_atoms)
                    {
                        if (a.monomer_ring_edges == null) a.monomer_ring_edges = new List<ring.ring_edge>();
                        a.monomer_ring_edges.Add(ring_monomer_edge);
                    }
                }

                foreach (var ring_monomer_node in ring_monomer_nodes)
                {
                    var node_atoms = atoms.Where(a => ring_monomer_node.Chain == a.chain_id && ring_monomer_node.Position == a.residue_index && ring_monomer_node.Residue1 == a.amino_acid).ToList();

                    foreach (var a in node_atoms)
                    {
                        if (a.monomer_ring_nodes == null) a.monomer_ring_nodes = new List<ring.ring_node>();
                         a.monomer_ring_nodes.Add(ring_monomer_node);
                    }
                }
            }
        }

        public static void load_dssp(string pdb_id, List<Atom> atoms)
        {
            pdb_id = Path.GetFileNameWithoutExtension(pdb_id);

            var multimer_dssp_list_file = $@"C:\betastrands_dataset\ss_dssp\{pdb_id}.dssp";

            var multimer_dssp_list = File.Exists(multimer_dssp_list_file) ? Dssp.Load(multimer_dssp_list_file).ToList() : new List<Dssp.DsspRecord>();

            var chain_ids = atoms.Select(a => a.chain_id).Distinct().ToList();


            var monomer_dssp_list = chain_ids.SelectMany(chain_id =>
            {
                var monomer_dssp_list_file = $@"C:\betastrands_dataset\ss_dssp\{pdb_id}{chain_id}.dssp";
                if (File.Exists(monomer_dssp_list_file))
                {
                    return Dssp.Load(monomer_dssp_list_file).ToList();
                }
                else
                {
                    return new List<Dssp.DsspRecord>();
                }
            }).ToList();

            //const char default_dssp_code = ' ';

            foreach (var atom in atoms)
            {
                {
                    var multimer_dssp = multimer_dssp_list.FirstOrDefault(a => (a.Chain == atom.chain_id) && (a.PdbResidueSequenceIndex.Trim() == atom.residue_index.ToString().Trim()) && (a.iCode == atom.i_code));

                    atom.multimer_dssp = multimer_dssp?.SecondaryStructure ?? Dssp.DsspRecord.dssp_default_secondary_structure;
                    atom.multimer_dssp3 = secondary_structure_state_reduction(atom.multimer_dssp);
                }
                {
                    var monomer_dssp = monomer_dssp_list.FirstOrDefault(a => (a.Chain == atom.chain_id) && (a.PdbResidueSequenceIndex.Trim() == atom.residue_index.ToString().Trim()) && (a.iCode == atom.i_code));

                    atom.monomer_dssp = monomer_dssp?.SecondaryStructure ?? Dssp.DsspRecord.dssp_default_secondary_structure;
                    atom.monomer_dssp3 = secondary_structure_state_reduction(atom.monomer_dssp);
                }
            }
        }


        public static void load_stride(string pdb_id, List<Atom> atoms)
        {
            pdb_id = Path.GetFileNameWithoutExtension(pdb_id);

            var stride_multimer_list_file = $@"C:\betastrands_dataset\ss_stride\{pdb_id}.stride";

            var stride_multimer_list = File.Exists(stride_multimer_list_file) ? StrideFile.Load(stride_multimer_list_file).Where(a => a.GetType() == typeof(StrideFile.Stride_DetailedSecondaryStructureAssignments)).Select(a => (StrideFile.Stride_DetailedSecondaryStructureAssignments) a).ToList() : new List<StrideFile.Stride_DetailedSecondaryStructureAssignments>();

            var chain_ids = atoms.Select(a => a.chain_id).Distinct().ToList();

            var stride_monomer_list = chain_ids.SelectMany(chain_id =>
            {
                var stride_monomer_list_file = $@"C:\betastrands_dataset\ss_stride\{pdb_id}{chain_id}.stride";

                if (File.Exists(stride_monomer_list_file))
                {
                    return StrideFile.Load(stride_monomer_list_file);
                }
                else
                {
                    return new List<StrideFile.Stride_Record>();
                }
            }).Where(a => a.GetType() == typeof(StrideFile.Stride_DetailedSecondaryStructureAssignments)).Select(a => (StrideFile.Stride_DetailedSecondaryStructureAssignments) a).ToList();

            foreach (var atom in atoms)
            {
                {
                    var multimer_stride = stride_multimer_list.FirstOrDefault(a => (a.ProteinChainIdentifier.Trim() == atom.chain_id.ToString().Trim()) && (a.PdbResidueNumber.Trim() == atom.residue_index.ToString().Trim()));

                    atom.multimer_stride = multimer_stride?.OneLetterSecondaryStructureCode?[0] ?? StrideFile.Stride_DetailedSecondaryStructureAssignments.stride_default_secondary_structure;
                    atom.multimer_stride3 = secondary_structure_state_reduction(atom.multimer_stride);
                }
                {
                    var monomer_stride = stride_monomer_list.FirstOrDefault(a => (a.ProteinChainIdentifier.Trim() == atom.chain_id.ToString().Trim()) && (a.PdbResidueNumber.Trim() == atom.residue_index.ToString().Trim()));

                    atom.monomer_stride = monomer_stride?.OneLetterSecondaryStructureCode?[0] ?? StrideFile.Stride_DetailedSecondaryStructureAssignments.stride_default_secondary_structure;
                    atom.monomer_stride3 = secondary_structure_state_reduction(atom.monomer_stride);
                }
            }
        }

        public static void load_rsa(string pdb_id, List<Atom> atoms)
        {
            pdb_id = Path.GetFileNameWithoutExtension(pdb_id);

            //atoms.Select(a => a.PdbId).ToList();

            var rsa_folder = @"C:\betastrands_dataset\sasa\";

            if (!Directory.Exists(rsa_folder)) return;

            var rsa_file_wildcard = $"{pdb_id}*.?-rsa-sasa";

            var rsa_files = Directory.GetFiles(rsa_folder, rsa_file_wildcard).ToList();
            var rsa_list = RSA.load(rsa_files);

            foreach (var atom in atoms)
            {
                atom.RSA_L = rsa_list.FirstOrDefault(a => (string.Equals(a.pdb_id, atom.pdb_id, StringComparison.InvariantCultureIgnoreCase)) && (a.chain_id == atom.chain_id) && (a.algo == 'L') && (a.res_num == atom.residue_index) && (a.amino_acid == atom.amino_acid));
                atom.RSA_S = rsa_list.FirstOrDefault(a => (string.Equals(a.pdb_id, atom.pdb_id, StringComparison.InvariantCultureIgnoreCase)) && (a.chain_id == atom.chain_id) && (a.algo == 'S') && (a.res_num == atom.residue_index) && (a.amino_acid == atom.amino_acid));
            }
        }

        public static double Distance3D(double x0, double y0, double z0, double x1, double y1, double z1)
        {
            return (double) Math.Sqrt((double) (((x1 - x0) * (x1 - x0)) + ((y1 - y0) * (y1 - y0)) + ((z1 - z0) * (z1 - z0))));
        }

        public static double measure_atomic_linear_distance(List<Atom> atoms)
        {
            atoms = select_amino_acid_master_atoms(null, atoms);

            var atom1 = atoms.First();
            var atom2 = atoms.Last();

            var distance = Distance3D(atom1.X, atom1.Y, atom1.Z, atom2.X, atom2.Y, atom2.Z);

            return distance;
        }

        public static double measure_atomic_curve_distance(List<Atom> atoms)
        {
            // atoms must be in correct order
            // measures the distance in angstrom of the curve
            atoms = select_amino_acid_master_atoms(null, atoms);

            var curve_distances = new double[atoms.Count - 1];
            for (var index = 0; index < atoms.Count - 1; index++)
            {
                var atom1 = atoms[index];
                var atom2 = atoms[index + 1];

                curve_distances[index] = Distance3D(atom1.X, atom1.Y, atom1.Z, atom2.X, atom2.Y, atom2.Z);
            }

            return curve_distances.Sum();
        }

        public static (double displacement, double distance_of_curve, double tortuosity1) measure_tortuosity1(List<Atom> atoms)
        {
            if (atoms == null || atoms.Count <= 1) return (0, 0, 0);


            atoms = select_amino_acid_master_atoms(null, atoms);

            var distance_of_curve = measure_atomic_curve_distance(atoms);

            //var distance_of_curve = distance_of_curve_stats.sum;

            var displacement = measure_atomic_linear_distance(atoms);

            var tortuosity1 = distance_of_curve / displacement;

            return (displacement, distance_of_curve, tortuosity1);
        }

        public static List<T> list_subsection_from_obj_to_obj<T>(List<T> list, T first_list_object, T last_list_object)
        {
            var index1 = list.IndexOf(first_list_object);
            var index2 = list.IndexOf(last_list_object);

            return list.Skip(index1 < index2 ? index1 : index2).Take(Math.Abs(index2 - index1) + 1).ToList();
        }

        public static (descriptive_stats displacement_stat_values, descriptive_stats curve_stat_values, descriptive_stats tortuosity_stat_values) measure_tortuosity2(List<Atom> atoms)
        {
            if (atoms == null || atoms.Count <= 1) return (null, null, null);

            atoms = select_amino_acid_master_atoms(null, atoms);

            var pairs = atoms.SelectMany((a, x) => atoms.Where((b, y) => x < y).Select((b, y) => (atom1: a, atom2: b)).ToList()).ToList();

            var distances = pairs.Select(a => (curve_distance: measure_atomic_curve_distance(list_subsection_from_obj_to_obj(atoms, a.atom1, a.atom2)), displacement_distance: Distance3D(a.atom1.X, a.atom1.Y, a.atom1.Z, a.atom2.X, a.atom2.Y, a.atom2.Z))).ToArray();

            // linear_distance / linear_list / linear_stat_values --> list of distances between all pairs of atoms, then averaged
            // curve_distance / 
            var displacement_distance_list = distances.Select(a => a.displacement_distance).ToArray(); // list of LineOfSight distances between all points (CA atoms)
            //var curve_stat_list = distances.Select(a => a.curve_distance).ToArray(); // list of all curve distances between all pairs of points (CA atoms) 
            var curve_distance_list = distances.Select(a => a.curve_distance).ToArray(); //curve_stat_list.Select(a => a.sum).ToArray();

            var tortuosity_list = distances.Select(a => a.curve_distance / a.displacement_distance).ToArray();

            descriptive_stats displacement_stat_values = descriptive_stats.get_stat_values(displacement_distance_list, nameof(displacement_stat_values)); // "linear_stat_values");
            descriptive_stats curve_stat_values = descriptive_stats.get_stat_values(curve_distance_list, nameof(curve_stat_values)); // "curve_stat_values");
            descriptive_stats tortuosity_stat_values = descriptive_stats.get_stat_values(tortuosity_list, nameof(tortuosity_stat_values)); //"tortuosity_stat_values");

            return (displacement_stat_values, curve_stat_values, tortuosity_stat_values);
        }

        public static void find_intramolecular_contacts(string pdb_id, List<Atom> atoms, double? max_dist = null)
        {
            const bool intermolecular = false;
            const bool intramolecular = true;
            find_contacts(pdb_id, atoms, intermolecular: intermolecular, intramolecular: intramolecular, max_dist: max_dist);
        }

        public static void find_intermolecular_contacts(string pdb_id, List<Atom> atoms, double? max_dist = null)
        {
            const bool intermolecular = true;
            const bool intramolecular = false;
            find_contacts(pdb_id, atoms, intermolecular: intermolecular, intramolecular: intramolecular, max_dist: max_dist);
        }

        private static void find_contacts(string pdb_id, List<Atom> atoms, bool intermolecular, bool intramolecular, double? max_dist = null)
        {
            var contact_element_types = new char[] {'C', 'N', 'O'};

            if (pdb_id != null) atoms = atoms.Where(a => string.Equals(a.pdb_id, pdb_id, StringComparison.InvariantCultureIgnoreCase)).ToList();

            if (intermolecular) atoms.ForEach(a => a.contact_map_intermolecular = new List<(Atom atom, double distance)>());
            if (intramolecular) atoms.ForEach(a => a.contact_map_intramolecular = new List<(Atom atom, double distance)>());

            if (contact_element_types != null && contact_element_types.Length > 0) atoms = atoms.Where(a => contact_element_types.Contains(a.element)).ToList();

            const bool skip_same_amino_acid_contacts = true;

            var pairs = atoms.SelectMany((a, x) => atoms.Where((b, y) => x < y).Select((b, y) => (atom1: a, atom2: b)).Where(b => (intramolecular && b.atom1.chain_id == b.atom2.chain_id && (!skip_same_amino_acid_contacts || b.atom1.residue_index != b.atom2.residue_index)) || (intermolecular && b.atom1.chain_id != b.atom2.chain_id)).ToList()).ToList();

            var pair_distances = pairs.Select(a => (atom1: a.atom1, atom2: a.atom2, distance: Distance3D(a.atom1.X, a.atom1.Y, a.atom1.Z, a.atom2.X, a.atom2.Y, a.atom2.Z))).ToList();

            pair_distances.ForEach(a =>
            {
                if (max_dist != null && a.distance > max_dist) return;


                if (a.atom1.chain_id == a.atom2.chain_id)
                {
                    a.atom1.contact_map_intramolecular.Add((a.atom2, a.distance));
                    a.atom2.contact_map_intramolecular.Add((a.atom1, a.distance));
                }

                if (a.atom1.chain_id != a.atom2.chain_id)
                {
                    a.atom1.contact_map_intermolecular.Add((a.atom2, a.distance));
                    a.atom2.contact_map_intermolecular.Add((a.atom1, a.distance));
                }
            });
        }


        public static List<(Atom atom1, Atom atom2, double distance)> get_master_contacts(string pdb_id, List<Atom> atoms)
        {
            var master_atoms = select_amino_acid_master_atoms(pdb_id, atoms);

            var pairs = master_atoms.SelectMany((a, x) => master_atoms.Where((b, y) => x < y).Select((b, y) => (atom1: a, atom2: b)).ToList()).ToList();

            var pair_distances = pairs.Select(a => (atom1: a.atom1, atom2: a.atom2, distance: Distance3D(a.atom1.X, a.atom1.Y, a.atom1.Z, a.atom2.X, a.atom2.Y, a.atom2.Z))).ToList();

            return pair_distances;

        }

        public static subsequence_classification_data get_intramolecular_protein_1d(subsequence_classification_data subsequence_data)
        {
            var result = new subsequence_classification_data()
            {
                pdb_id = subsequence_data.pdb_id,
                chain_id = subsequence_data.chain_id,
                class_id = subsequence_data.class_id,

                dimer_type = subsequence_data.dimer_type,
                parallelism = subsequence_data.parallelism,
                symmetry_mode = subsequence_data.symmetry_mode,

                aa_subsequence = null,
                res_ids = null,

                dssp_monomer_subsequence = null,
                dssp_multimer_subsequence = null,

                stride_monomer_subsequence = null,
                stride_multimer_subsequence = null,

                subsequence_atoms = subsequence_data.pdb_chain_atoms,
                subsequence_master_atoms = subsequence_data.pdb_chain_master_atoms,

                pdb_chain_atoms = subsequence_data.pdb_chain_atoms,
                pdb_chain_master_atoms = subsequence_data.pdb_chain_master_atoms,

                parent = subsequence_data,

                foldx_energy_differences = null
            };

            
            
            var nh_sequence_list = Atom.amino_acid_sequence(subsequence_data.pdb_id, result.subsequence_master_atoms);
            var nh_dssp_multimer_list = Atom.ss_sequence(subsequence_data.pdb_id, result.subsequence_master_atoms, structure_oligomisation.multimer, feature_calcs.ss_types.DSSP);
            var nh_dssp_monomer_list = Atom.ss_sequence(subsequence_data.pdb_id, result.subsequence_master_atoms, structure_oligomisation.monomer, feature_calcs.ss_types.DSSP);
            var nh_stride_multimer_list = Atom.ss_sequence(subsequence_data.pdb_id, result.subsequence_master_atoms, structure_oligomisation.multimer, feature_calcs.ss_types.STRIDE);
            var nh_stride_monomer_list = Atom.ss_sequence(subsequence_data.pdb_id, result.subsequence_master_atoms, structure_oligomisation.monomer, feature_calcs.ss_types.STRIDE);

            result.aa_subsequence = (nh_sequence_list.Count > 0) ? nh_sequence_list.First().aa_sequence : null;
            result.res_ids = result.subsequence_master_atoms.Select(a => (a.residue_index, a.i_code, a.amino_acid)).ToList();
            result.dssp_monomer_subsequence = (nh_dssp_monomer_list.Count > 0) ? nh_dssp_monomer_list.First().ss_sequence : null;
            result.dssp_multimer_subsequence = (nh_dssp_multimer_list.Count > 0) ? nh_dssp_multimer_list.First().ss_sequence : null;
            result.stride_monomer_subsequence = (nh_stride_monomer_list.Count > 0) ? nh_stride_monomer_list.First().ss_sequence : null;
            result.stride_multimer_subsequence = (nh_stride_multimer_list.Count > 0) ? nh_stride_multimer_list.First().ss_sequence : null;
            
            return result;
        }

        public static subsequence_classification_data get_intramolecular_neighbourhood_1d(subsequence_classification_data subsequence_data, int neighbourhood_flanking_size = 6)
        {
            if (neighbourhood_flanking_size % 2 != 0 || neighbourhood_flanking_size % 3 != 0)
            {
                throw new Exception("Number must be divisible by both 2 and 3");
            }

            var result = new subsequence_classification_data()
            {
                pdb_id = subsequence_data.pdb_id,
                chain_id = subsequence_data.chain_id,
                class_id = subsequence_data.class_id,

                dimer_type = subsequence_data.dimer_type,
                parallelism = subsequence_data.parallelism,
                symmetry_mode = subsequence_data.symmetry_mode,

                aa_subsequence = null,
                res_ids = null,

                dssp_monomer_subsequence = null,
                dssp_multimer_subsequence = null,

                stride_monomer_subsequence = null,
                stride_multimer_subsequence = null,

                subsequence_atoms = null,
                subsequence_master_atoms = null,

                pdb_chain_atoms = subsequence_data.pdb_chain_atoms,
                pdb_chain_master_atoms = subsequence_data.pdb_chain_master_atoms,

                parent = subsequence_data,

                foldx_energy_differences = null
            };


            // returns the neighbours of a primary amino acid sequence
            var subsequence_atoms_array_indexes = subsequence_data.subsequence_master_atoms.Select(a => subsequence_data.pdb_chain_master_atoms.IndexOf(a)).Where(a => a != -1).ToList();
            var start_array_index = subsequence_atoms_array_indexes.Min(); // first res_id in interface or subseq
            var end_array_index = subsequence_atoms_array_indexes.Max(); // last res_id in interface or subseq

            var all_array_indexes = subsequence_data.subsequence_master_atoms.Select(a => a.residue_index).ToList();
            var total_indexes_before_start = all_array_indexes.Count(a => a < start_array_index);
            var total_indexes_after_end = all_array_indexes.Count(a => a > end_array_index);

            var neighbourhood_size_before_start = neighbourhood_flanking_size;
            var neighbourhood_size_after_end = neighbourhood_flanking_size;

            if (total_indexes_before_start < neighbourhood_flanking_size && total_indexes_after_end > neighbourhood_flanking_size)
            {
                neighbourhood_size_after_end = neighbourhood_size_after_end + (neighbourhood_flanking_size - total_indexes_before_start);
            }

            if (total_indexes_after_end < neighbourhood_flanking_size && total_indexes_before_start > neighbourhood_flanking_size)
            {
                neighbourhood_size_before_start = neighbourhood_size_before_start + (neighbourhood_flanking_size - total_indexes_after_end);
            }

            result.subsequence_master_atoms = subsequence_data.pdb_chain_master_atoms.Where((a, i) => (i >= start_array_index - neighbourhood_size_before_start && i < start_array_index) || (i > end_array_index && i <= end_array_index + neighbourhood_size_after_end)).Except(subsequence_data.subsequence_atoms).ToList();


            result.subsequence_atoms = result.subsequence_master_atoms.SelectMany(a => a.amino_acid_atoms).Distinct().ToList();

            result.res_ids = result.subsequence_master_atoms.Select(a => (a.residue_index, a.i_code, a.amino_acid)).ToList();

            var nh_sequence_list = Atom.amino_acid_sequence(subsequence_data.pdb_id, result.subsequence_master_atoms);
            var nh_dssp_multimer_list = Atom.ss_sequence(subsequence_data.pdb_id, result.subsequence_master_atoms, structure_oligomisation.multimer, feature_calcs.ss_types.DSSP);
            var nh_dssp_monomer_list = Atom.ss_sequence(subsequence_data.pdb_id, result.subsequence_master_atoms, structure_oligomisation.monomer, feature_calcs.ss_types.DSSP);
            var nh_stride_multimer_list = Atom.ss_sequence(subsequence_data.pdb_id, result.subsequence_master_atoms, structure_oligomisation.multimer, feature_calcs.ss_types.STRIDE);
            var nh_stride_monomer_list = Atom.ss_sequence(subsequence_data.pdb_id, result.subsequence_master_atoms, structure_oligomisation.monomer, feature_calcs.ss_types.STRIDE);

            result.aa_subsequence = (nh_sequence_list.Count > 0) ? nh_sequence_list.First().aa_sequence : null;
            result.dssp_multimer_subsequence = (nh_dssp_multimer_list.Count > 0) ? nh_dssp_multimer_list.First().ss_sequence : null;
            result.dssp_monomer_subsequence = (nh_dssp_monomer_list.Count > 0) ? nh_dssp_monomer_list.First().ss_sequence : null;
            result.stride_multimer_subsequence = (nh_stride_multimer_list.Count > 0) ? nh_stride_multimer_list.First().ss_sequence : null;
            result.stride_monomer_subsequence = (nh_stride_monomer_list.Count > 0) ? nh_stride_monomer_list.First().ss_sequence : null;

            
            return result;
        }

        public static subsequence_classification_data get_intramolecular_protein_3d(subsequence_classification_data subsequence_data)
        {
            var result = new subsequence_classification_data()
            {
                pdb_id = subsequence_data.pdb_id,
                chain_id = subsequence_data.chain_id,
                class_id = subsequence_data.class_id,

                dimer_type = subsequence_data.dimer_type,
                parallelism = subsequence_data.parallelism,
                symmetry_mode = subsequence_data.symmetry_mode,

                aa_subsequence = null,
                res_ids = null,

                dssp_monomer_subsequence = null,
                dssp_multimer_subsequence = null,

                stride_monomer_subsequence = null,
                stride_multimer_subsequence = null,

                subsequence_atoms = subsequence_data.pdb_chain_atoms,
                subsequence_master_atoms = subsequence_data.pdb_chain_master_atoms,

                pdb_chain_atoms = subsequence_data.pdb_chain_atoms,
                pdb_chain_master_atoms = subsequence_data.pdb_chain_master_atoms,

                parent = subsequence_data,

                foldx_energy_differences = null
            };            

            var protein_sequence_list = Atom.amino_acid_sequence(subsequence_data.pdb_id, result.subsequence_master_atoms);
            var protein_dssp_multimer_list = Atom.ss_sequence(subsequence_data.pdb_id, result.subsequence_master_atoms, structure_oligomisation.multimer, feature_calcs.ss_types.DSSP);
            var protein_dssp_monomer_list = Atom.ss_sequence(subsequence_data.pdb_id, result.subsequence_master_atoms, structure_oligomisation.monomer, feature_calcs.ss_types.DSSP);
            var protein_stride_multimer_list = Atom.ss_sequence(subsequence_data.pdb_id, result.subsequence_master_atoms, structure_oligomisation.multimer, feature_calcs.ss_types.STRIDE);
            var protein_stride_monomer_list = Atom.ss_sequence(subsequence_data.pdb_id, result.subsequence_master_atoms, structure_oligomisation.monomer, feature_calcs.ss_types.STRIDE);

            result.aa_subsequence = (protein_sequence_list.Count > 0) ? protein_sequence_list.First().aa_sequence : null;
            result.res_ids = result.subsequence_master_atoms.Select(a => (a.residue_index, a.i_code, a.amino_acid)).ToList();
            result.dssp_monomer_subsequence = (protein_dssp_monomer_list.Count > 0) ? protein_dssp_monomer_list.First().ss_sequence : null;
            result.dssp_multimer_subsequence = (protein_dssp_multimer_list.Count > 0) ? protein_dssp_multimer_list.First().ss_sequence : null;
            result.stride_monomer_subsequence = (protein_stride_monomer_list.Count > 0) ? protein_stride_monomer_list.First().ss_sequence : null;
            result.stride_multimer_subsequence = (protein_stride_multimer_list.Count > 0) ? protein_stride_multimer_list.First().ss_sequence : null;

            var fx = foldx_caller.calc_energy_differences(result.pdb_id, result.chain_id, result.res_ids, false, subsequence_classification_data.protein_data_sources.protein_3d);
            result.foldx_energy_differences = fx;

            // todo: check fx isn't empty for neighbourhood / protein

            return result;
        }

        public static subsequence_classification_data get_intramolecular_neighbourhood_3d(subsequence_classification_data subsequence_data, double max_dist = (double) 5.0)
        {
            if (max_dist > 8.0)
            {
                throw new Exception("max atomic distance too long for contat");
            }

            var result = new subsequence_classification_data()
            {
                pdb_id = subsequence_data.pdb_id,
                chain_id = subsequence_data.chain_id,
                class_id = subsequence_data.class_id,

                dimer_type = subsequence_data.dimer_type,
                parallelism = subsequence_data.parallelism,
                symmetry_mode = subsequence_data.symmetry_mode,

                aa_subsequence = null,
                res_ids = null,

                dssp_monomer_subsequence = null,
                dssp_multimer_subsequence = null,

                stride_monomer_subsequence = null,
                stride_multimer_subsequence = null,

                subsequence_atoms = null,
                subsequence_master_atoms = null,

                pdb_chain_atoms = subsequence_data.pdb_chain_atoms,
                pdb_chain_master_atoms = subsequence_data.pdb_chain_master_atoms,

                parent = subsequence_data,
            };

            // returns the 3d contact neighbours of a sequence
            var all_3d_contacts = subsequence_data.subsequence_atoms.SelectMany(a => a.amino_acid_atoms.SelectMany(b => b.contact_map_intramolecular).ToList()).Distinct().OrderBy(a => a.distance).ToList();

            var neighbourhood_3d_contacts = all_3d_contacts.Where(a => a.distance <= max_dist).Select(a => a.atom.amino_acid_master_atom).Distinct().ToList(); //.Take(10).ToList();

            var min_3d_nh_length = 3;

            if (neighbourhood_3d_contacts.Count < min_3d_nh_length)
            {
                neighbourhood_3d_contacts = all_3d_contacts.Select(a => a.atom.amino_acid_master_atom).Distinct().Take(min_3d_nh_length).ToList();
            }

            neighbourhood_3d_contacts = neighbourhood_3d_contacts.Select(a => a.amino_acid_master_atom).Distinct().ToList();

            result.subsequence_atoms = neighbourhood_3d_contacts;

            result.subsequence_master_atoms = Atom.select_amino_acid_master_atoms(result.pdb_id, result.subsequence_atoms).Except(subsequence_data.subsequence_atoms).ToList();
            result.res_ids = result.subsequence_master_atoms.Select(a => (a.residue_index, a.i_code, a.amino_acid)).ToList();

            var nh_sequence_list        = Atom.amino_acid_sequence(subsequence_data.pdb_id, result.subsequence_master_atoms);
            var nh_dssp_multimer_list   = Atom.ss_sequence(subsequence_data.pdb_id, result.subsequence_master_atoms, structure_oligomisation.multimer, feature_calcs.ss_types.DSSP);
            var nh_dssp_monomer_list    = Atom.ss_sequence(subsequence_data.pdb_id, result.subsequence_master_atoms, structure_oligomisation.monomer, feature_calcs.ss_types.DSSP);
            var nh_stride_multimer_list = Atom.ss_sequence(subsequence_data.pdb_id, result.subsequence_master_atoms, structure_oligomisation.multimer, feature_calcs.ss_types.STRIDE);
            var nh_stride_monomer_list  = Atom.ss_sequence(subsequence_data.pdb_id, result.subsequence_master_atoms, structure_oligomisation.monomer, feature_calcs.ss_types.STRIDE);

            result.aa_subsequence = (nh_sequence_list.Count > 0) ? nh_sequence_list.First().aa_sequence : null;
            result.dssp_multimer_subsequence = (nh_dssp_multimer_list.Count > 0) ? nh_dssp_multimer_list.First().ss_sequence : null;
            result.dssp_monomer_subsequence = (nh_dssp_monomer_list.Count > 0) ? nh_dssp_monomer_list.First().ss_sequence : null;
            result.stride_multimer_subsequence = (nh_stride_multimer_list.Count > 0) ? nh_stride_multimer_list.First().ss_sequence : null;
            result.stride_monomer_subsequence = (nh_stride_monomer_list.Count > 0) ? nh_stride_monomer_list.First().ss_sequence : null;

            var fx = foldx_caller.calc_energy_differences(result.pdb_id, result.chain_id, result.res_ids, false, subsequence_classification_data.protein_data_sources.neighbourhood_3d);
            result.foldx_energy_differences = fx;

            // todo: check fx isn't empty for neighbourhood / protein
            return result;
        }
    }
}