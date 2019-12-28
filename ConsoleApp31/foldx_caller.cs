using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Diagnostics.Contracts;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Threading.Tasks;

namespace Coil_DHC_Dataset_Maker
{
    public static class foldx_caller
    {
        public const string foldx_folder = @"C:\betastrands_dataset\foldx\";
        public const string pdb_folder = @"C:\betastrands_dataset\foldx\pdb\";



        //public class energy_differences
        //{
        //    public List<(bool is_fragment, bool is_repaired, double[] energy)> ala_scan;
        //    public List<(bool is_fragment, bool is_repaired, double[] energy)> build_model_position_scan;
        //    public List<(bool is_fragment, bool is_repaired, double[] energy)> build_model_subsequence_mutant;
        //    public List<(bool is_fragment, bool is_repaired, double[] energy)> position_scan;
        //    //public List<(bool is_fragment, bool is_repaired, double[] energy)> stability;
        //}

        public class energy_differences
        {
            //public (string cmd_line, string wait_filename, List<foldx_ala_scanning_result> data) foldx_ala_scanning_result_protein;

            public (string cmd_line, string wait_filename, List<foldx_ala_scanning_result> data) foldx_ala_scanning_result_subsequence;

            public (string cmd_line, string wait_filename, List<foldx_position_scanning_result> data) foldx_position_scanning_result_subsequence;

            public (string cmd_line, string wait_filename, List<foldx_energy_terms_ps> data) foldx_buildmodel_position_scan_result_subsequence;

            public (string cmd_line, string wait_filename, List<foldx_energy_terms_sm> data) foldx_buildmodel_subsequence_mutant_result_subsequence;
        }


        private static object file_write_lock = new object();

        public static energy_differences calc_energy_differences(string pdb_id, char chain_id, List<(int residue_index, char i_code, char amino_acid)> res_ids, bool run, subsequence_classification_data.protein_data_sources source)//int nh_first_res_id, int nh_last_res_id)
        {
#if DEBUG
            if (Program.verbose_debug) Program.WriteLine($"{nameof(calc_energy_differences)}(string pdb_id, char chain_id, List<(int residue_index, char i_code, char amino_acid)> res_ids, bool run);");
#endif

            //var nh_first_res_id = res_ids.Min(a => a.residue_index);
            //var nh_last_res_id = res_ids.Max(a => a.residue_index);

            var result = new energy_differences();
            // fragments should be first residue to last residue (even non-interacting residues), since we want to mutate only some of them

            // returns energy differences upon mutation of a neighbourhood
            // can be viewed as: difference over whole monomer, difference over fragment only, and also, repaired/unrepaired, but we only do repaired.

            //var pdb_folder = @"c:\betastrands_dataset\pdb\";
            //var foldx_folder = @"C:\betastrands_dataset\foldx\";
            pdb_id = Path.GetFileNameWithoutExtension(pdb_id).Substring(0, 4);

            // make monomer from dimer
            var monomer_file = Atom.extract_split_pdb_chains(pdb_id, chain_id).First();

            // repair
            var monomer_file_repair = foldx_repair_pdb(Path.GetFileNameWithoutExtension(monomer_file), run);

            var repair_res_ids = File.ReadAllLines(monomer_file_repair).Where(a => a.StartsWith("ATOM")).Select(a => int.Parse(a.Substring(22, 4))).Distinct().OrderBy(a => a).ToList();

            // filter res ids to remove res ides which were in the original pdb structure file but removed by the foldx repair
            res_ids = res_ids.Where(a => repair_res_ids.Contains(a.residue_index)).ToList();

            // calculate energy differences before/after mutation

            //result.foldx_ala_scanning_result_protein = foldx_caller.load_foldx_ala_scanning(monomer_file_repair, chain_id, null, run);
            result.foldx_ala_scanning_result_subsequence = foldx_caller.load_foldx_ala_scanning(monomer_file_repair, chain_id, res_ids, run);
            result.foldx_position_scanning_result_subsequence             = foldx_caller.load_foldx_position_scanning((monomer_file_repair, chain_id, res_ids), run);
            result.foldx_buildmodel_position_scan_result_subsequence      = foldx_caller.load_foldx_buildmodel_position_scan((monomer_file_repair, chain_id, res_ids), run);
            result.foldx_buildmodel_subsequence_mutant_result_subsequence = foldx_caller.load_foldx_buildmodel_subsequence_mutant((monomer_file_repair, chain_id, res_ids), run);

            lock (file_write_lock)
            {
                File.AppendAllLines(Path.Combine($@"{foldx_folder}", $"foldx_calc_ala_scanning_{source.ToString()}.bat.skip"), new[]
                {
                    $@"if not exist ""{result.foldx_ala_scanning_result_subsequence.wait_filename}"" {result.foldx_ala_scanning_result_subsequence.cmd_line}"
                });
                File.AppendAllLines(Path.Combine($@"{foldx_folder}", $"foldx_calc_position_scanning_{source.ToString()}.bat"), new[]
                {
                    $@"if not exist ""{result.foldx_position_scanning_result_subsequence.wait_filename}"" {result.foldx_position_scanning_result_subsequence.cmd_line}"
                });
                File.AppendAllLines(Path.Combine($@"{foldx_folder}", $"foldx_calc_buildmodel_position_scan_{source.ToString()}.bat"), new[]
                {
                    $@"if not exist ""{result.foldx_buildmodel_position_scan_result_subsequence.wait_filename}"" {result.foldx_buildmodel_position_scan_result_subsequence.cmd_line}"
                });
                File.AppendAllLines(Path.Combine($@"{foldx_folder}", $"foldx_calc_buildmodel_subsequence_mutant_{source.ToString()}.bat"), new[]
                {
                    $@"if not exist ""{result.foldx_buildmodel_subsequence_mutant_result_subsequence.wait_filename}"" {result.foldx_buildmodel_subsequence_mutant_result_subsequence.cmd_line}"
                });
            }
        

            return result;
        }



        public class foldx_energy_terms
        {
            public string pdb_id;
            public char chain_id;
            public List<(int residue_index, char i_code, char amino_acid)> res_ids;

            public int line_index;

            //public string pdb_filename;
            //public string wait_filename;
            //public bool repaired;
            public string Pdb;
            public double SD;
            public double total_energy;
            public double Backbone_Hbond;
            public double Sidechain_Hbond;
            public double Van_der_Waals;
            public double Electrostatics;
            public double Solvation_Polar;
            public double Solvation_Hydrophobic;
            public double Van_der_Waals_clashes;
            public double entropy_sidechain;
            public double entropy_mainchain;
            public double sloop_entropy;
            public double mloop_entropy;
            public double cis_bond;
            public double torsional_clash;
            public double backbone_clash;
            public double helix_dipole;
            public double water_bridge;
            public double disulfide;
            public double electrostatic_kon;
            public double partial_covalent_bonds;
            public double energy_Ionisation;
            public double Entropy_Complex;

            public (string name, double value)[] properties() {
                return new (string name, double value)[]
                    {
                        (nameof(total_energy),              total_energy),
                        (nameof(Backbone_Hbond),            Backbone_Hbond),
                        (nameof(Sidechain_Hbond),           Sidechain_Hbond),
                        (nameof(Van_der_Waals),             Van_der_Waals),
                        (nameof(Electrostatics),            Electrostatics),
                        (nameof(Solvation_Polar),           Solvation_Polar),
                        (nameof(Solvation_Hydrophobic),     Solvation_Hydrophobic),
                        (nameof(Van_der_Waals_clashes),     Van_der_Waals_clashes),
                        (nameof(entropy_sidechain),         entropy_sidechain),
                        (nameof(entropy_mainchain),         entropy_mainchain),
                        (nameof(sloop_entropy),             sloop_entropy),
                        (nameof(mloop_entropy),             mloop_entropy),
                        (nameof(cis_bond),                  cis_bond),
                        (nameof(torsional_clash),           torsional_clash),
                        (nameof(backbone_clash),            backbone_clash),
                        (nameof(helix_dipole),              helix_dipole),
                        (nameof(water_bridge),              water_bridge),
                        (nameof(disulfide),                 disulfide),
                        (nameof(electrostatic_kon),         electrostatic_kon),
                        (nameof(partial_covalent_bonds),    partial_covalent_bonds),
                        (nameof(energy_Ionisation),         energy_Ionisation),
                        (nameof(Entropy_Complex),           Entropy_Complex),
                    };
            }
        }


        public static readonly
            List<(int index, string full_name, string foldx_aa_code3, char foldx_aa_code1, string standard_aa_code3, char standard_aa_code1, bool is_mutable, string residue_type)>
            foldx_residues =
                new
                    List<(int index, string full_name, string foldx_aa_code3, char foldx_aa_code1, string standard_aa_code3, char standard_aa_code1, bool is_mutable, string residue_type)>()
                    { // add MSE/MSA ?? the special MET substitute

                        (00, "glycine", "GLY", 'G', "GLY", 'G', true, "standard aa"),
                        (01, "alanine", "ALA", 'A', "ALA", 'A', true, "standard aa"),
                        (02, "leucine", "LEU", 'L', "LEU", 'L', true, "standard aa"),
                        (03, "valine", "VAL", 'V', "VAL", 'V', true, "standard aa"),
                        (04, "isoleucine", "ILE", 'I', "ILE", 'I', true, "standard aa"),
                        (05, "proline", "PRO", 'P', "PRO", 'P', true, "standard aa"),
                        (06, "arginine", "ARG", 'R', "ARG", 'R', true, "standard aa"),
                        (07, "threonine", "THR", 'T', "THR", 'T', true, "standard aa"),
                        (08, "serine", "SER", 'S', "SER", 'S', true, "standard aa"),
                        (09, "cysteine", "CYS", 'C', "CYS", 'C', true, "standard aa"),
                        (10, "methionine", "MET", 'M', "MET", 'M', true, "methionine"),
                        (11, "lysine", "LYS", 'K', "LYS", 'K', true, "standard aa"),
                        (12, "glutamic", "GLU", 'E', "GLU", 'E', true, "standard aa"),
                        (13, "glutamine", "GLN", 'Q', "GLN", 'Q', true, "standard aa"),
                        (14, "aspartic", "ASP", 'D', "ASP", 'D', true, "standard aa"),
                        (15, "asparagine", "ASN", 'N', "ASN", 'N', true, "standard aa"),
                        (16, "tryptophane", "TRP", 'W', "TRP", 'W', true, "standard aa"),
                        (17, "tyrosine", "TYR", 'Y', "TYR", 'Y', true, "standard aa"),
                        (18, "phenylalanine", "PHE", 'F', "PHE", 'F', true, "standard aa"),
                        (19, "histidine", "HIS", 'H', "HIS", 'H', true, "standard aa"),
                        (20, "phoshoporylated threonine", "PTR", 'y', "THR", 'T', true, "phoshoporylated threonine"),
                        (21, "phosphorylated tyrosine", "TPO", 'p',  "TYR", 'Y', true, "phosphorylated tyrosine"),
                        (22, "phosphorylated serine", "SEP", 's', "SER", 'S', true, "phosphorylated serine"),
                        (23, "hydroxiproline", "HYP", 'h', "PRO", 'P', true, "hydroxiproline"),
                        (24, "sulfotyrosine", "TYS", 'z', "TYR", 'Y', true, "sulfotyrosine"),
                        (25, "monomethylated lysine", "MLZ", 'k', "LYS", 'K', true, "monomethylated lysine"),
                        (26, "dimethylated lysine", "MLY", 'm', "LYS", 'K', true, "dimethylated lysine"),
                        (27, "trimethylated lysine", "M3L", 'l', "LYS", 'K', true, "trimethylated lysine"),
                        (28, "charged ND1 histidine", "H1S", 'o', "HIS", 'H', true, "charged ND1 histidine"),
                        (29, "charged NE2 histidine", "H2S", 'e', "HIS", 'H', true, "charged NE2 histidine"),
                        (30, "neutral histidine", "H3S", 'f', "HIS", 'H', true, "neutral histidine"),
                        (31, "adenosine", "A", 'a', "A", 'a', true, "adenosine"),
                        (32, "guanosine", "G", 'g', "G", 'g', true, "guanosine"),
                        (33, "cytosine", "C", 'c', "C", 'c', true, "cytosine"),
                        (34, "thymidine", "T", 't',"T", 't', true, "thymidine"),
                        (35, "6-methylated adenosine", "6MA", 'b', "6MA", 'b', true, "6-methylated adenosine"),
                        (36, "5-methylated cytosine", "5CM", 'd', "5CM", 'd', true, "5-methylated cytosine"),
                        (37, "adenosine triphosphate", "ATP", ' ', "ATP", ' ', false, "adenosine triphosphate"),
                        (38, "adenosine diphosphate", "ADP", ' ', "ADP", ' ', false, "adenosine diphosphate"),
                        (39, "guanosine triphosphate", "GTP", ' ', "GTP", ' ', false, "guanosine triphosphate"),
                        (40, "guanosine diphosphate", "GDP", ' ',"GDP", ' ', false, "guanosine diphosphate"),
                        (41, "dioctadecylglycerol-3-phosphatidyl-choline", "LIP", ' ',"LIP", ' ', false, "dioctadecylglycerol-3-phosphatidyl-choline"),
                        (42, "calcium", "CA", ' ',"CA", ' ', false, "calcium"),
                        (43, "magnesium", "MG", ' ',"MG", ' ', false, "magnesium"),
                        (44, "manganese", "MN", ' ',"MN", ' ', false, "manganese"),
                        (45, "sodium", "NA", ' ',"NA", ' ', false, "sodium"),
                        (46, "zinc", "ZN", ' ',"ZN", ' ', false, "zinc"),
                        (47, "iron", "FE", ' ',"FE", ' ', false, "iron"),
                        (48, "copper", "CU", ' ',"CU", ' ', false, "copper"),
                        (49, "cobalt", "CO", ' ',"CO", ' ', false, "cobalt"),
                        (50, "potasium", "K", ' ',"K", ' ', false, "potasium"),
                        (51, "water", "HOH", ' ',"HOH", ' ', false, "water"),
                    };

        public static readonly List<(int index, string full_name, string foldx_aa_code3, char foldx_aa_code1, string standard_aa_code3, char standard_aa_code1, bool is_mutable, string residue_type)>
            foldx_residues_aa_mutable = foldx_residues.Where(a => (a.index >= 0) && (a.index <= 30)).ToList();

        private static readonly List<string> call_foldx_lock = new List<string>();

        public static string[] read_all_lines_until_success(string wait_file, int max_tries = int.MaxValue, int delay_ms = 10)
        {
#if DEBUG
            if (Program.verbose_debug) Program.WriteLine($"{nameof(read_all_lines_until_success)}(string wait_file, int max_tries = int.MaxValue, int delay_ms = 10);");
#endif

            var file_accessed = false;
            var access_attempt_count = 0;
            string[] data = null;

            while (!file_accessed)
            {
                if (access_attempt_count >= max_tries) throw new IOException($"File could not be accessed \"{wait_file}\"");

                access_attempt_count++;

                try
                {

                    if (File.Exists(wait_file))// && new FileInfo(wait_file).Length > 0)
                    {
                        var f_len1 = new FileInfo(wait_file).Length;

                        data = File.ReadAllLines(wait_file);

                        var f_len2 = new FileInfo(wait_file).Length;

                        if (f_len1 != f_len2)
                        {
                            data = null;

                            Task.Delay(delay_ms).Wait();

                        }
                        else
                        {
                            file_accessed = true;

                        }
                    }
                    else
                    {
                        Task.Delay(delay_ms).Wait();

                    }

                }
                catch (IOException e)
                {
                    lock (Program._console_lock)
                    {
                        Console.WriteLine(e.Message);
                    }

                    Task.Delay(delay_ms).Wait();
                }
            }

            return data;
        }

        public static (string cmd_line, string[] data) call_foldx(string pdb_file_id, string foldx_command, string foldx_args, string wait_filename, string lock_code_local, bool run, string pdb_folder = foldx_caller.pdb_folder)
        {
#if DEBUG
            if (Program.verbose_debug) Program.WriteLine($@"call_foldx(pdb_file_id={pdb_file_id}, foldx_command={foldx_command}, foldx_args={foldx_args}, wait_filename={wait_filename}, lock_code_local={lock_code_local}");
#endif

            var foldx_exe = $@"{foldx_folder}foldx.exe";
            var pdb_file = $"{pdb_folder}{Path.GetFileNameWithoutExtension(pdb_file_id)}.pdb";


            var lock_code = string.Join("_", new string[] { pdb_file_id, foldx_command, foldx_args, wait_filename, lock_code_local });


            var wait_for_complete = false;

            lock (call_foldx_lock)
            {
                if (call_foldx_lock.Contains(lock_code))
                {
                    wait_for_complete = true;
                }
                else
                {
                    call_foldx_lock.Add(lock_code);
                }
            }

            while (wait_for_complete)
            {
                Task.Delay(10).Wait();

                lock (call_foldx_lock)
                {
                    wait_for_complete = call_foldx_lock.Contains(lock_code);

                    if (!wait_for_complete)
                    {
                        call_foldx_lock.Add(lock_code);
                    }
                }
            }

            string[] data = null;


            var foldx_args2 = foldx_args;
            if (!foldx_args2.Contains("--command=")) { foldx_args2 += $" --command={foldx_command}"; }
            if (!foldx_args2.Contains("--pdb=")) { foldx_args2 += $" --pdb={Path.GetFileName(pdb_file)}"; }
            if (!foldx_args2.Contains("--pdb-dir=")) { foldx_args2 += $" --pdb-dir=\"{Path.GetDirectoryName(pdb_file) ?? "."}\""; }
            if (!foldx_args2.Contains("--out-pdb=") && !string.Equals(foldx_command, "RepairPDB", StringComparison.InvariantCultureIgnoreCase)) { foldx_args2 += $" --out-pdb=false"; }

            var start = new ProcessStartInfo { FileName = foldx_exe, WorkingDirectory = Path.GetDirectoryName(foldx_exe) ?? "", Arguments = foldx_args2, UseShellExecute = false, CreateNoWindow = false, RedirectStandardOutput = true, RedirectStandardError = true };
            var cmd_line = $@"""{start.FileName}"" {start.Arguments}";

            if (string.IsNullOrWhiteSpace(wait_filename) || !File.Exists(wait_filename) || new FileInfo(wait_filename).Length <= 0)
            {

                if (run)
                {

                    lock (Program._console_lock)
                    {
                        Console.WriteLine($"{nameof(call_foldx)}: run: \"" + start.FileName + "\" " + start.Arguments);
                    }

                    using (var process = Process.Start(start))
                    {
                        if (process == null) throw new Exception($"{nameof(call_foldx)}: {nameof(process)} is null");


                        using (var reader = process.StandardOutput)
                        {
                            var stdout = reader.ReadToEnd();
                            stdout = stdout.Replace("\r\n", "\r\nstdout: ");
                            if (!string.IsNullOrWhiteSpace(stdout))
                            {
                                lock (Program._console_lock)
                                {
                                    Console.WriteLine($"{nameof(stdout)}: {stdout}");
                                }
                            }

                            var stderr = process.StandardError.ReadToEnd();
                            stderr = stderr.Replace("\r\n", "\r\nstderr: ");

                            if (!string.IsNullOrWhiteSpace(stderr))
                            {
                                lock (Program._console_lock)
                                {
                                    Console.WriteLine($"{nameof(call_foldx)}: {stderr}");
                                }
                            }
                        }

                        process.WaitForExit();
                    }

                    if (run)
                    {
                        data = read_all_lines_until_success(wait_filename);
                    }
                    else
                    {
                        data = File.Exists(wait_filename) && new FileInfo(wait_filename).Length > 0 ? File.ReadAllLines(wait_filename) : new string[0];
                    }
                }
            }
            else
            {
                if (!string.IsNullOrWhiteSpace(wait_filename))
                {
                    if (run)
                    {
                        data = read_all_lines_until_success(wait_filename);
                    }
                    else
                    {
                        data = File.Exists(wait_filename) && new FileInfo(wait_filename).Length > 0 ? File.ReadAllLines(wait_filename) : new string[0];
                    }
                }
            }


            lock (call_foldx_lock)
            {
                call_foldx_lock.Remove(lock_code);
            }

            return (cmd_line, data);
        }

        public static string foldx_repair_pdb(string pdb_id, bool run, string pdb_folder = foldx_caller.pdb_folder, string repair_folder = @"C:\betastrands_dataset\foldx\pdb\")
        {
#if DEBUG
            if (Program.verbose_debug) Program.WriteLine($"foldx_repair_pdb(pdb_id = \"{pdb_id}\")");
#endif
            var pdb_file =        Path.Combine($@"{pdb_folder}", $@"{Path.GetFileNameWithoutExtension(pdb_id)}.pdb");
            var pdb_file_repair = Path.Combine($@"{repair_folder}", $@"{Path.GetFileNameWithoutExtension(pdb_id)}_Repair.pdb");

            if (File.Exists(pdb_file_repair) && new FileInfo(pdb_file_repair).Length > 0) return pdb_file_repair;

            var foldx_cmd = "RepairPDB";
            var foldx_args = $"--output-dir=\"{Path.GetDirectoryName(pdb_file) ?? "."}\"";
            var call_foldx_result = call_foldx(pdb_id, foldx_cmd, foldx_args, pdb_file_repair, pdb_id, run, pdb_folder);

            return pdb_file_repair;
        }


        public class foldx_energy_terms_ps : foldx_energy_terms
        {
            public (char original_amino_acid1, char chain_id, int residue_index, char mutant_foldx_amino_acid1, string mutant_foldx_amino_acid3, char mutant_standard_amino_acid1, string mutant_standard_amino_acid3) mutation_positions_data;
        }

        public static (string cmd_line, string wait_filename, List<foldx_energy_terms_ps> data) load_foldx_buildmodel_position_scan((string pdb_id, char chain_id, List<(int residue_index, char i_code, char amino_acid)> res_ids) interface_residues, bool run)
        {
            var (pdb_id, chain_id, res_ids) = interface_residues;

#if DEBUG
            if (Program.verbose_debug) Program.WriteLine($"load_foldx_buildmodel_position_scan(pdb_id = \"{pdb_id}\", chain_id = \"{chain_id}\", res_ids = \"{res_ids}\")");
#endif

            if (res_ids == null || res_ids.Count == 0) return ("", "", null);

            pdb_id = Path.GetFileNameWithoutExtension(pdb_id);
            var lock_code = $"bm_ps_{pdb_id}{chain_id}{string.Join("_", res_ids)}";

            var foldx_mutation_positions_data = new List<string>();

            var mutation_positions_data = new List<(
                char original_amino_acid, 
                char chain_id,
                int residue_index, 
                char mutant_foldx_amino_acid1, 
                string mutant_foldx_amino_acid3, 
                char mutant_standard_amino_acid1,
                string mutant_standard_amino_acid3)>();

            foreach (var master_atom in res_ids)
            {
                foreach (var mutant_amino_acid in foldx_residues_aa_mutable)//.Select(b => b.foldx_aa_code1).ToList())
                {

                    var mp_data = (
                        original_amino_acid:master_atom.amino_acid, 
                        chain_id,
                        master_atom.residue_index, 
                        mutant_foldx_amino_acid1:mutant_amino_acid.foldx_aa_code1, 
                        mutant_foldx_amino_acid3:mutant_amino_acid.foldx_aa_code3, 
                        mutant_standard_amino_acid1:mutant_amino_acid.standard_aa_code1, 
                        mutant_standard_amino_acid3:mutant_amino_acid.standard_aa_code3);

                    mutation_positions_data.Add(mp_data);

                    var mp = $"{mp_data.original_amino_acid}{mp_data.chain_id}{mp_data.residue_index}{mp_data.mutant_foldx_amino_acid1};";
                    foldx_mutation_positions_data.Add(mp);
                }
            }

            var first_amino_acid = $"{res_ids.First().amino_acid}{res_ids.First().residue_index}";
            var last_amino_acid = $"{res_ids.Last().amino_acid}{res_ids.Last().residue_index}";
            var reside_index_sum = res_ids.Select(a => a.residue_index).Sum();

            

            var mutant_list_file = Path.Combine(foldx_folder, "bm_ps", $@"individual_list_bm_ps_{pdb_id}{chain_id}_{first_amino_acid}_{last_amino_acid}_{reside_index_sum}.txt");

            lock (file_write_lock)
            {
                Directory.CreateDirectory(Path.GetDirectoryName(mutant_list_file));
                File.WriteAllLines(mutant_list_file, foldx_mutation_positions_data);
            }

            var foldx_cmd = $"BuildModel";
            var foldx_output_pdb = false;
            var foldx_number_of_runs = 1;
            var output_file_tag = $"bm_ps_{pdb_id}_{first_amino_acid}_{last_amino_acid}_{reside_index_sum}";
            var wait_filename = Path.Combine(foldx_folder, "bm_ps", $"Dif_{output_file_tag}_{pdb_id}.fxout");
            var foldx_args = $"--mutant-file={mutant_list_file} --numberOfRuns={foldx_number_of_runs} --out-pdb={foldx_output_pdb} --output-file={output_file_tag}";

            var foldx_result_1 = call_foldx(pdb_id, foldx_cmd, foldx_args, wait_filename, lock_code, run);
            var cmd_line = foldx_result_1.cmd_line;
            var foldx_result = foldx_result_1.data;

            if (foldx_result == null || foldx_result.Length == 0)
            {
                return (cmd_line, wait_filename, new List<foldx_energy_terms_ps>());
            }

            var sd = false;
            var marker_index = foldx_result.ToList().FindIndex(a =>a.StartsWith("Pdb\ttotal energy"));

            if (marker_index < 0)
            {
                marker_index = foldx_result.ToList().FindIndex(a =>a.StartsWith("Pdb\tSD\ttotal energy"));

                if (marker_index > -1)
                {
                    sd = true;
                }
            }

            var foldx_result2 = foldx_result.Skip(marker_index + 1).ToList();

            var results = foldx_result2.Select((a, i) =>
            {
                var b = a.Split();

                var j = 0;

                var c = new foldx_energy_terms_ps()
                {
                    line_index = i,
                    pdb_id = pdb_id,
                    chain_id = chain_id,
                    res_ids = res_ids,

                    mutation_positions_data = mutation_positions_data[i],

                    Pdb = b[j++],
                    SD = sd?double.Parse(b[j++]):0,
                    total_energy = double.Parse(b[j++]),
                    Backbone_Hbond = double.Parse(b[j++]),
                    Sidechain_Hbond = double.Parse(b[j++]),
                    Van_der_Waals = double.Parse(b[j++]),
                    Electrostatics = double.Parse(b[j++]),
                    Solvation_Polar = double.Parse(b[j++]),
                    Solvation_Hydrophobic = double.Parse(b[j++]),
                    Van_der_Waals_clashes = double.Parse(b[j++]),
                    entropy_sidechain = double.Parse(b[j++]),
                    entropy_mainchain = double.Parse(b[j++]),
                    sloop_entropy = double.Parse(b[j++]),
                    mloop_entropy = double.Parse(b[j++]),
                    cis_bond = double.Parse(b[j++]),
                    torsional_clash = double.Parse(b[j++]),
                    backbone_clash = double.Parse(b[j++]),
                    helix_dipole = double.Parse(b[j++]),
                    water_bridge = double.Parse(b[j++]),
                    disulfide = double.Parse(b[j++]),
                    electrostatic_kon = double.Parse(b[j++]),
                    partial_covalent_bonds = double.Parse(b[j++]),
                    energy_Ionisation = double.Parse(b[j++]),
                    Entropy_Complex = double.Parse(b[j++])

                };

                return c;
            }).ToList();

            if (results.Count != foldx_mutation_positions_data.Count)
            {
                throw new Exception("Missing foldx model results");
            }

            return (cmd_line, wait_filename, results);
        }


        public class foldx_energy_terms_sm : foldx_energy_terms
        {
            public List<(char original_amino_acid1, char chain_id, int residue_index, char mutant_foldx_amino_acid1, string mutant_foldx_amino_acid3, char mutant_standard_amino_acid1, string mutant_standard_amino_acid3)> mutation_positions_data;
        }

        public static (string cmd_line, string wait_filename, List<foldx_energy_terms_sm> data) load_foldx_buildmodel_subsequence_mutant((string pdb_id, char chain_id, List<(int residue_index, char i_code, char amino_acid)> res_ids) interface_residues, bool run)
        {
            var (pdb_id, chain_id, res_ids) = interface_residues;

            if (res_ids == null || res_ids.Count == 0) return ("", "", null);

#if DEBUG
            if (Program.verbose_debug) Program.WriteLine($"load_foldx_buildmodel_subsequence_mutant(pdb_id = \"{pdb_id}\", chain_id = \"{chain_id}\", res_ids = \"{res_ids}\")");
#endif
            // FoldX --command=BuildModel --pdb=BM.pdb --mutant-file=individual_list.txt
            // individual_list.txt = e.g. FA39L,FB39L; (res_aa, chain, res_num, mutant_aa)


            pdb_id = Path.GetFileNameWithoutExtension(pdb_id);
            var lock_code = $"bm_if_subs_{pdb_id}{chain_id}{string.Join("_", res_ids)}";

            var mutation_positions_data = new List<List<(char original_amino_acid1, char chain_id, int residue_index, char mutant_foldx_amino_acid1, string mutant_foldx_amino_acid3, char mutant_standard_amino_acid1, string mutant_standard_amino_acid3)>>();
            var foldx_mutation_positions_data = new List<string>();

            foreach (var mutant_amino_acid in foldx_residues_aa_mutable)
            {
                var mp_data = res_ids.Select(a => (original_amino_acid1:a.amino_acid, chain_id:chain_id, residue_index:a.residue_index, mutant_foldx_amino_acid1:mutant_amino_acid.foldx_aa_code1, mutant_foldx_amino_acid3:mutant_amino_acid.foldx_aa_code3, mutant_standard_amino_acid1:mutant_amino_acid.standard_aa_code1, mutant_standard_amino_acid3:mutant_amino_acid.standard_aa_code3)).ToList();
                mutation_positions_data.Add(mp_data);

                var mp = $"{string.Join(",", mp_data.Select(a => $"{a.original_amino_acid1}{a.chain_id}{a.residue_index}{a.mutant_foldx_amino_acid1}").ToList())};";

                foldx_mutation_positions_data.Add(mp);
            }

            var first_amino_acid = $"{res_ids.First().amino_acid}{res_ids.First().residue_index}";
            var last_amino_acid = $"{res_ids.Last().amino_acid}{res_ids.Last().residue_index}";
            var reside_index_sum = res_ids.Select(a => a.residue_index).Sum();

            var mutant_list_file = Path.Combine( foldx_folder, "bm_if_subs", $@"individual_list_bm_if_subs_{pdb_id}{chain_id}_{first_amino_acid}_{last_amino_acid}_{reside_index_sum}.txt");

            lock (file_write_lock)
            {
                Directory.CreateDirectory(Path.GetDirectoryName(mutant_list_file));
                File.WriteAllLines(mutant_list_file, foldx_mutation_positions_data);
            }


            var foldx_cmd = $"BuildModel";
            var foldx_output_pdb = false;
            var foldx_number_of_runs = 1;
            var output_file_tag = $"bm_if_subs_{pdb_id}_{first_amino_acid}_{last_amino_acid}_{reside_index_sum}";
            var wait_filename = Path.Combine(foldx_folder, "bm_if_subs", $@"Dif_{output_file_tag}_{pdb_id}.fxout");

            var foldx_args = $"--mutant-file={mutant_list_file} --numberOfRuns={foldx_number_of_runs} --out-pdb={foldx_output_pdb} --output-file={output_file_tag}";

            var foldx_result_1 = call_foldx(pdb_id, foldx_cmd, foldx_args, wait_filename, lock_code, run);
            var cmd_line = foldx_result_1.cmd_line;

            var foldx_result = foldx_result_1.data;

            if (foldx_result == null || foldx_result.Length <= 0)
            {
                return (cmd_line, wait_filename, new List<foldx_energy_terms_sm>());
            }

            var sd = false;
            var marker_index = foldx_result.ToList().FindIndex(a => a.StartsWith("Pdb\ttotal energy"));

            if (marker_index < 0)
            {
                marker_index = foldx_result.ToList().FindIndex(a => a.StartsWith("Pdb\tSD\ttotal energy"));

                if (marker_index > -1)
                {
                    sd = true;
                }
            }

            var foldx_result2 = foldx_result.Skip(marker_index + 1).ToList();

            var results = foldx_result2.Select((a, i) =>
            {
                var b = a.Split();

                var j = 0;

                var foldx_energy_terms = new foldx_energy_terms_sm
                {
                    mutation_positions_data = mutation_positions_data[i],
                    line_index = i,
                    pdb_id = pdb_id,
                    chain_id = chain_id,
                    res_ids = res_ids,

                    Pdb = b[j++],
                    SD = sd ? double.Parse(b[j++]) : 0,
                    total_energy = double.Parse(b[j++]),
                    Backbone_Hbond = double.Parse(b[j++]),
                    Sidechain_Hbond = double.Parse(b[j++]),
                    Van_der_Waals = double.Parse(b[j++]),
                    Electrostatics = double.Parse(b[j++]),
                    Solvation_Polar = double.Parse(b[j++]),
                    Solvation_Hydrophobic = double.Parse(b[j++]),
                    Van_der_Waals_clashes = double.Parse(b[j++]),
                    entropy_sidechain = double.Parse(b[j++]),
                    entropy_mainchain = double.Parse(b[j++]),
                    sloop_entropy = double.Parse(b[j++]),
                    mloop_entropy = double.Parse(b[j++]),
                    cis_bond = double.Parse(b[j++]),
                    torsional_clash = double.Parse(b[j++]),
                    backbone_clash = double.Parse(b[j++]),
                    helix_dipole = double.Parse(b[j++]),
                    water_bridge = double.Parse(b[j++]),
                    disulfide = double.Parse(b[j++]),
                    electrostatic_kon = double.Parse(b[j++]),
                    partial_covalent_bonds = double.Parse(b[j++]),
                    energy_Ionisation = double.Parse(b[j++]),
                    Entropy_Complex = double.Parse(b[j++])
                };

                return foldx_energy_terms;

            }).ToList();

            if (results.Count != foldx_mutation_positions_data.Count)
            {
                throw new Exception("Missing foldx model results");
            }


            return (cmd_line, wait_filename, results);
        }

        public class foldx_ala_scanning_result
        {
            // LYS 86 to ALA energy change is 0.0306484
            // THR 87 to ALA energy change is -0.727468
            // GLU 88 to ALA energy change is -0.923282

            public string pdb_id;
            public char chain_id;
            public int residue_index;

            //public string res_name;
            //public string scan_res_name;

            public char original_foldx_amino_acid_1;
            public string original_foldx_amino_acid_3;
            public char original_standard_amino_acid_1;
            public string original_standard_amino_acid_3;

            public char mutant_foldx_amino_acid_1;
            public string mutant_foldx_amino_acid_3;
            public char mutant_standard_amino_acid_1;
            public string mutant_standard_amino_acid_3;

            public double ddg;
        }

        public class foldx_position_scanning_result
        {
            // FoldX --command=PositionScan --pdb-dir=c:\betastrands_dataset\pdb\ --pdb=1AJYA.pdb --positions MA30d

            // META30M 0
            // META30G 0.0440333
            // META30A 0.000395327
            // META30L -0.063261

            public string pdb_id;
            public char chain_id;
            public int residue_index;

            public char original_foldx_amino_acid_1;
            public string original_foldx_amino_acid_3;
            public char original_standard_amino_acid_1;
            public string original_standard_amino_acid_3;

            public char mutant_foldx_amino_acid_1;
            public string mutant_foldx_amino_acid_3;
            public char mutant_standard_amino_acid_1;
            public string mutant_standard_amino_acid_3;

            public double ddg;
        }

        public static string fix_non_standard_naming(string res_name, bool fix)
        {
            if (fix)
            {
                if (res_name.Length == 3)
                {
                    var ri = foldx_residues.FindIndex(a => a.foldx_aa_code3 == res_name);
                    if (ri > -1)
                    {
                        var r = foldx_residues[ri];

                        if (r.foldx_aa_code3 != r.standard_aa_code3)
                        {
                            return r.standard_aa_code3;
                        }

                        return res_name;
                    }
                }

                if (res_name.Length == 1)
                {
                    var ri = foldx_residues.FindIndex(a => a.foldx_aa_code1 == res_name[0]);
                    if (ri > -1)
                    {
                        var r = foldx_residues[ri];

                        if (r.foldx_aa_code1 != r.standard_aa_code1)
                        {
                            return r.standard_aa_code1.ToString(CultureInfo.InvariantCulture);
                        }

                        return res_name;
                    }
                }
            }

            return res_name;
        }

        public static (string cmd_line, string wait_filename, List<foldx_ala_scanning_result> data) load_foldx_ala_scanning(string pdb_id, char chain_id, List<(int residue_index, char i_code, char amino_acid)> res_ids, bool run)
        {
            //    var (pdb_id, chain_id, res_ids) = interface_residues;
#if DEBUG
            //if (Program.verbose_debug) Program.WriteLine($"load_foldx_ala_scanning(pdb_id = \"{pdb_id}\", chain_id = \"{chain_id}\", res_ids = \"{res_ids}\")");
#endif
            pdb_id = Path.GetFileNameWithoutExtension(pdb_id);

            //if (res_ids == null || res_ids.Count == 0) return ("", "", null);

            var wait_filename = Path.Combine(foldx_folder,"ala_scan",$@"{pdb_id}_AS.fxout");

            string[] file_data = null;



            //var cmd_line = "";


            //if (!File.Exists(wait_filename))
            //{

            //var first_amino_acid = $"{res_ids.First().amino_acid}{res_ids.First().residue_index}";
            //var last_amino_acid = $"{res_ids.Last().amino_acid}{res_ids.Last().residue_index}";
            //var reside_index_sum = res_ids.Select(a => a.residue_index).Sum();

            var lock_code = $"{pdb_id}{chain_id}{(res_ids != null ? string.Join("_", res_ids) : "")}";

            var foldx_cmd = $"AlaScan";
            var foldx_args = ""; //$"--output-file={file_tag}";

            var call_foldx_result_1 = call_foldx(pdb_id, foldx_cmd, foldx_args, wait_filename, lock_code, run);
            var cmd_line = call_foldx_result_1.cmd_line;

            var call_foldx_result = call_foldx_result_1.data;

            file_data = call_foldx_result;

            //}
            //else
            //{
            //    file_data = File.ReadAllLines(wait_filename);
            //}

            if (file_data == null || file_data.Length == 0)
            {
                return (cmd_line, wait_filename, new List<foldx_ala_scanning_result>());
            }

            //if (file_data != null && file_data.Length > 0)
            // {
            var split_lines = file_data.Select(b => b.Split()).ToList();

            //var fix_nsn = false;

            var results = split_lines.Select(c => 
                new foldx_ala_scanning_result()
                {
                    pdb_id = pdb_id,
                    chain_id = chain_id,
                    residue_index = int.Parse(c[1]),

                    //res_name = fix_non_standard_naming(c[0], fix_nsn),

                    original_foldx_amino_acid_1 = foldx_residues_aa_mutable.First(a =>string.Equals(a.foldx_aa_code3, c[0], StringComparison.InvariantCultureIgnoreCase)).foldx_aa_code1,
                    original_standard_amino_acid_1 = foldx_residues_aa_mutable.First(a =>string.Equals(a.foldx_aa_code3, c[0], StringComparison.InvariantCultureIgnoreCase)).standard_aa_code1,
                    original_foldx_amino_acid_3 = foldx_residues_aa_mutable.First(a =>string.Equals(a.foldx_aa_code3, c[0], StringComparison.InvariantCultureIgnoreCase)).foldx_aa_code3,
                    original_standard_amino_acid_3 = foldx_residues_aa_mutable.First(a =>string.Equals(a.foldx_aa_code3, c[0], StringComparison.InvariantCultureIgnoreCase)).standard_aa_code3,

                    mutant_foldx_amino_acid_1 = foldx_residues_aa_mutable.First(a =>string.Equals(a.foldx_aa_code3, c[3], StringComparison.InvariantCultureIgnoreCase)).foldx_aa_code1,
                    mutant_standard_amino_acid_1 = foldx_residues_aa_mutable.First(a =>string.Equals(a.foldx_aa_code3, c[3], StringComparison.InvariantCultureIgnoreCase)).standard_aa_code1,
                    mutant_foldx_amino_acid_3 = foldx_residues_aa_mutable.First(a =>string.Equals(a.foldx_aa_code3, c[3], StringComparison.InvariantCultureIgnoreCase)).foldx_aa_code3,
                    mutant_standard_amino_acid_3 = foldx_residues_aa_mutable.First(a =>string.Equals(a.foldx_aa_code3, c[3], StringComparison.InvariantCultureIgnoreCase)).standard_aa_code3,


                    ddg = double.Parse(c[7], NumberStyles.Float, CultureInfo.InvariantCulture)
                }).ToList();

            

            results = results.Where(a => pdb_id == a.pdb_id && chain_id == a.chain_id).ToList();

            if (res_ids != null && res_ids.Count > 0)
            {
                results = results.Where(a => res_ids.Any(b => b.residue_index == a.residue_index)).ToList();
            }
            

            return (cmd_line, wait_filename, results);
        }



        public static (string cmd_line, string wait_filename, List<foldx_position_scanning_result> data) load_foldx_position_scanning((string pdb_id, char chain_id, List<(int residue_index, char i_code, char amino_acid)> res_ids) interface_residues, bool run)
        {
            var (pdb_id, chain_id, res_ids) = interface_residues;

#if DEBUG
            if (Program.verbose_debug) Program.WriteLine($"load_foldx_position_scanning(pdb_id = \"{pdb_id}\", chain_id = \"{chain_id}\", res_ids = \"{res_ids}\")");
#endif

            if (res_ids == null || res_ids.Count == 0) return ("", "", null);

            pdb_id = Path.GetFileNameWithoutExtension(pdb_id);

            var lock_code = $"{pdb_id}{chain_id}{string.Join("_", res_ids)}";

            var first_amino_acid = $"{res_ids.First().amino_acid}{res_ids.First().residue_index}";
            var last_amino_acid = $"{res_ids.Last().amino_acid}{res_ids.Last().residue_index}";
            var reside_index_sum = res_ids.Select(a => a.residue_index).Sum();

            var mutation_code = 'd';
            var mutation_position_ids = res_ids.Select(a => $"{a.amino_acid}{chain_id}{a.residue_index}{mutation_code}").ToList();
            var mutation_positions = $"{string.Join(",", mutation_position_ids)}";
            var file_tag = $"{pdb_id}_{mutation_position_ids.First()}_{mutation_position_ids.Last()}_{reside_index_sum}";
            var wait_filename = Path.Combine(foldx_folder, "ps", $@"PS_{file_tag}_scanning_output.txt");
            var foldx_cmd = $"PositionScan";
            var foldx_args = $"--positions={mutation_positions} --output-file={file_tag}";

            var call_foldx_result_1 = call_foldx(pdb_id, foldx_cmd, foldx_args, wait_filename, lock_code, run);
            var call_foldx_result = call_foldx_result_1.data;

            List<foldx_position_scanning_result> results = null;

            //var fix_nsn = false;

            if (call_foldx_result != null && call_foldx_result.Length > 0)
            {
                var split_lines = call_foldx_result.Select(b => b.Split()).ToList();

                results = split_lines.Select(line =>
                {
                    
                    return new foldx_position_scanning_result()
                    {
                        pdb_id = pdb_id,
                        chain_id = chain_id,
                        

                        original_foldx_amino_acid_1 = foldx_residues_aa_mutable.First(a=> string.Equals(a.foldx_aa_code3, line[0].Substring(0, 3), StringComparison.InvariantCultureIgnoreCase)).foldx_aa_code1,
                        original_standard_amino_acid_1 = foldx_residues_aa_mutable.First(a => string.Equals(a.foldx_aa_code3, line[0].Substring(0, 3), StringComparison.InvariantCultureIgnoreCase)).standard_aa_code1,

                        original_foldx_amino_acid_3 = foldx_residues_aa_mutable.First(a => string.Equals(a.foldx_aa_code3, line[0].Substring(0, 3), StringComparison.InvariantCultureIgnoreCase)).foldx_aa_code3,
                        original_standard_amino_acid_3 = foldx_residues_aa_mutable.First(a => string.Equals(a.foldx_aa_code3, line[0].Substring(0, 3), StringComparison.InvariantCultureIgnoreCase)).standard_aa_code3,

                        mutant_foldx_amino_acid_1 = foldx_residues_aa_mutable.First(a => a.foldx_aa_code1.Equals(line[0][line[0].Length - 1])).foldx_aa_code1,
                        mutant_standard_amino_acid_1 = foldx_residues_aa_mutable.First(a => a.foldx_aa_code1.Equals(line[0][line[0].Length - 1])).standard_aa_code1,
                        mutant_foldx_amino_acid_3 = foldx_residues_aa_mutable.First(a => a.foldx_aa_code1.Equals(line[0][line[0].Length - 1])).foldx_aa_code3,
                        mutant_standard_amino_acid_3 = foldx_residues_aa_mutable.First(a => a.foldx_aa_code1.Equals(line[0][line[0].Length - 1])).standard_aa_code3,

                        residue_index = int.Parse(line[0].Substring(4, line[0].Length - 5)),
                        //res_mutant_aa = line[0][line[0].Length - 1],
                        ddg = double.Parse(line[1], NumberStyles.Float, CultureInfo.InvariantCulture)
                    };
                }).ToList();
            }

            return (call_foldx_result_1.cmd_line, wait_filename, results);
        }
    }
}
