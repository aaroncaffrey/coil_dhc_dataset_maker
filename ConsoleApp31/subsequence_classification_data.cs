using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Net;
using System.Net.Configuration;
using System.Net.Mime;
using System.Runtime.InteropServices;
using System.Security.Cryptography.X509Certificates;
using System.Security.Policy;
using System.Threading;
using System.Threading.Tasks;
using static Coil_DHC_Dataset_Maker.Program;
using Newtonsoft.Json;

namespace Coil_DHC_Dataset_Maker
{
    public partial class subsequence_classification_data
    {
        // xtodo: pssm feature(s) - mutatability

        // xtodo: foldx alanine scanning for energy change measure features

        // xtodo: monomeric sequence secondary structure prediction features

        // todo: uniprot composition feature

        // not_todo: alpha/beta protein percentage feature (alraedy included by dssp/stride feature)

        // subsequence properties
        public int class_id;
        public string dimer_type;
        public string parallelism;
        public string symmetry_mode;

        public string pdb_id;
        public char chain_id;
        public List<(int residue_index, char i_code, char amino_acid)> res_ids;
        //public List<(int res_id, char i_code, char res_aa)> neighbourhood_1d_res_ids;
        //public List<(int res_id, char i_code, char res_aa)> neighbourhood_3d_res_ids;

        public string aa_subsequence;
        //public string dssp_multimer_sequence;
        public string dssp_monomer_subsequence;
        //public string stride_multimer_sequence;
        public string stride_monomer_subsequence;

        public foldx_caller.energy_differences foldx_energy_differences;

        public string dssp_multimer_subsequence;
        public string stride_multimer_subsequence;
        public List<Atom> pdb_chain_atoms;
        public List<Atom> pdb_chain_master_atoms;
        public List<Atom> subsequence_atoms;
        public List<Atom> subsequence_master_atoms;

        //public descriptive_stats intramolecular_contact_count;
        //public descriptive_stats intramolecular_contact_distance;
        //public descriptive_stats foldx_monomer_ala_scan_ddg;
        //public descriptive_stats foldx_monomer_unrepaired_ala_scan_ddg;
        //public descriptive_stats foldx_monomer_repaired_ala_scan_ddg;
        //public double foldx_monomer_unrepaired_build_model_ddg;
        //public double foldx_monomer_repaired_build_model_ddg;

        public subsequence_classification_data parent = null;
        public subsequence_classification_data neighbourhood_1d = null;
        public subsequence_classification_data neighbourhood_3d = null;
        public subsequence_classification_data protein_1d = null;
        public subsequence_classification_data protein_3d = null;

        public static feature_info calculate_class_id_classification_data(subsequence_classification_data scd)
        {
#if DEBUG
            if (Program.verbose_debug) Program.WriteLine($"{nameof(calculate_class_id_classification_data)}(subsequence_classification_data scd);");
#endif


            var class_id_feature = new feature_info()
            {
                alphabet = "Overall",
                dimension = 0,
                category = nameof(scd.class_id),
                source = nameof(scd.class_id),
                group = nameof(scd.class_id),
                member = nameof(scd.class_id),
                perspective = nameof(scd.class_id),
                feature_value = scd.class_id
            };


            return class_id_feature;
        }

        public static List<feature_info> calculate_sable_classification_data_template = null;

        public static List<feature_info> calculate_sable_sequence_classification_data(subsequence_classification_data scd, List<Atom> subsequence_master_atoms, protein_data_sources source, int max_features)
        {
#if DEBUG
            if (Program.verbose_debug) Program.WriteLine($"{nameof(calculate_sable_sequence_classification_data)}(subsequence_classification_data scd, List<Atom> subsequence_master_atoms, protein_data_sources source, int max_features);");
#endif

            if (subsequence_master_atoms == null || subsequence_master_atoms.Count == 0)
            {
                if (calculate_sable_classification_data_template == null) throw new Exception();

                var template = calculate_sable_classification_data_template.Select(a => new feature_info(a)
                {
                    source = source.ToString(),
                    feature_value = 0
                }).ToList();

                return template;
            }

            var features = new List<feature_info>();

            //var pdb_sable_data = scd.pdb_chain_master_atoms.Select(a => a.sable_item).ToList();
            var subseq_sable_data = subsequence_master_atoms.Select(a => a.sable_item).ToList();

            var sable_data_list = new List<(string name, List<sable.sable_item> data)>();
            //sable_data_list.Add((nameof(pdb_sable_data), pdb_sable_data));
            sable_data_list.Add((nameof(subseq_sable_data), subseq_sable_data));


            foreach (var sable_data in sable_data_list)
            {
                var ds_entropy = descriptive_stats.get_stat_values(sable_data.data.Select(a => a.entropy_value).ToArray(), nameof(sable.sable_item.entropy_value));
                var ds_burial_abs = descriptive_stats.get_stat_values(sable_data.data.Select(a => a.absolute_burial_value).ToArray(), nameof(sable.sable_item.absolute_burial_value));
                var ds_burial_rel = descriptive_stats.get_stat_values(sable_data.data.Select(a => a.relative_burial_value).ToArray(), nameof(sable.sable_item.relative_burial_value));

                var feats = new List<feature_info>();

                var x0 = ds_entropy.encode().Select(a => new feature_info()
                {
                    alphabet = "Overall",
                    category = "sable",
                    dimension = 1,
                    source = source.ToString(),
                    @group = "sable_entropy",
                    member = a.member_id,
                    perspective = a.perspective_id,
                    feature_value = a.perspective_value
                }).ToList();

                var x1 = ds_burial_abs.encode().Select(a => new feature_info()
                {
                    alphabet = "Overall",
                    category = "sable",
                    dimension = 1,
                    source = source.ToString(),
                    @group = "sable_burial_abs",
                    member = a.member_id,
                    perspective = a.perspective_id,
                    feature_value = a.perspective_value
                }).ToList();

                var x2 = ds_burial_rel.encode().Select(a => new feature_info()
                {
                    alphabet = "Overall",
                    category = "sable",
                    dimension = 1,
                    source = source.ToString(),
                    @group = "sable_burial_rel",
                    member = a.member_id,
                    perspective = a.perspective_id,
                    feature_value = a.perspective_value
                }).ToList();

                feats.AddRange(x0);
                feats.AddRange(x1);
                feats.AddRange(x2);

                var x3 = feats.Select(a => new feature_info(a)
                {
                    @group = "sable_all"
                }).ToList();

                feats.AddRange(x3);

                features.AddRange(feats);
            }

            if (calculate_sable_classification_data_template == null)
            {
                calculate_sable_classification_data_template = features.Select(a => new feature_info(a) { source = "", feature_value = 0 }).ToList();
            }

            return features;
        }


        public static List<feature_info> calculate_mpsa_classification_data_template = null;
        public static List<feature_info> calculate_mpsa_classification_data(List<Atom> subsequence_master_atoms, protein_data_sources source, int max_features)
        {
#if DEBUG
            if (Program.verbose_debug) Program.WriteLine($"{nameof(calculate_mpsa_classification_data)}(List<Atom> subsequence_master_atoms, protein_data_sources source, int max_features);");
#endif

            if (subsequence_master_atoms == null || subsequence_master_atoms.Count == 0)
            {
                if (calculate_mpsa_classification_data_template == null) throw new Exception();

                var template = calculate_mpsa_classification_data_template.Select(a => new feature_info(a)
                {
                    source = source.ToString(),
                    feature_value = 0
                }).ToList();

                return template;
            }


            var complete_sequence = subsequence_master_atoms;



            var split_sequence = feature_calcs.split_sequence(complete_sequence, 3, 0, false);
            var sequences = new List<(string name, List<Atom> sequence)>();
            sequences.Add(("unsplit", complete_sequence));
            sequences.AddRange(split_sequence.Select(a => ("split", a)).ToList());

            var features = new List<feature_info>();


            for (var sq_index = 0; sq_index < sequences.Count; sq_index++)
            {
                var sq = sequences[sq_index];
                
                var features_dist_all = new List<feature_info>();
                var features_prob_all = new List<feature_info>();
                var features_dist_aa_all = new List<feature_info>();
                var features_prob_aa_all = new List<feature_info>();

                var master_indexes = sq.sequence.Select(a => a.master_index).Distinct().ToList();

                var mpsa_readers = sq.sequence.SelectMany(a => a.mpsa_entries.Select(b => b.mpsa_entry.reader).Distinct().ToList()).Distinct().ToList();
                //mpsa_readers = mpsa_readers.Where(a => a != null && a.mpsa_matrix != null && a.mpsa_matrix.Count > 0).ToList();
                mpsa_readers = mpsa_readers.Select(a => new mpsa_reader(a, master_indexes)).ToList();

                // add consensus mpsa_reader (modify mpsa_readers collection to change consensus source)
                mpsa_readers.Add(new mpsa_reader(mpsa_readers));


                var atoms_aa_seq = string.Join("", sq.sequence.Select(a => a.amino_acid).ToList());
                foreach (var reader in mpsa_readers)
                {
                    var format = reader.format;

                    var mpsa_aa_seq = string.Join("", reader.mpsa_matrix.Select(a => a.amino_acid).ToList());

                    if (atoms_aa_seq != mpsa_aa_seq) throw new Exception();

                    var ss_seq = string.Join("", reader.mpsa_matrix.Select(a => a.predicted_ss_code).ToList());

                    if (string.Equals(sq.name, "unsplit", StringComparison.InvariantCultureIgnoreCase))
                    {
                        // PseSSC data (secondary structure sequence pattern descriptors)
                        var mpsa_pse_aac_options = new pse_aac_options()
                        {
                            //composition = true, motifs = false, order_context = false, order_distance = true, position = true, split_composition = true

                            oaac = true,
                            oaac_binary = true,
                            motifs = true, //false
                            motifs_binary = true, //false
                            dipeptides = true, //false
                            dipeptides_binary = true, //false
                            //saac = true,
                            //saac_binary = true,
                            average_seq_position = true,
                            average_dipeptide_distance = true,
                        };

                        var pse_ssc = calculate_aa_or_ss_sequence_classification_data(source, "mpsa", $"mpsa_{sq.name}_{format}", ss_seq, feature_calcs.seq_type.secondary_structure_sequence, mpsa_pse_aac_options, max_features);

                        foreach (var a in pse_ssc.GroupBy(a => (a.alphabet, a.dimension, a.category, a.source, a.@group)).Select(feature_infos => feature_infos.ToList()).Where(a => a.Count <= max_features))
                        {
                            features.AddRange(a);
                        }
                    }

                    // Probability Data (probability of each amino acid being Helix, Strand, Coil or Turn)

                    var ss_overall_average = reader.ss_overall_average;
                    //var split_ss_overall_average = reader.split_ss_overall_average;

                    var merged_ss_overall_average = new List<(string group_suffix, string member_suffix, List<(char ss, double prob_value, double dist_value)> data)>();
                    merged_ss_overall_average.Add(( /*"unsplit"*/"", "", ss_overall_average));
                    //merged_ss_overall_average.AddRange(split_ss_overall_average.Select((a, i) => ("split", $"{i}", a)).ToList());

                    foreach (var item in merged_ss_overall_average)
                    {
                        foreach (var x in item.data)
                        {
                            var prob = new feature_info()
                            {
                                alphabet = "Overall",
                                dimension = 1,
                                category = "mpsa",
                                source = source.ToString(),
                                @group = $@"mpsa_{sq.name}_{format}_overall_probability_{item.group_suffix}",
                                member = $"{sq_index}_{x.ss}_{item.member_suffix}",
                                perspective = "default",
                                feature_value = x.prob_value,
                            };
                            features.Add(prob);

                            var dist = new feature_info()
                            {
                                alphabet = "Overall",
                                dimension = 1,
                                category = "mpsa",
                                source = source.ToString(),
                                @group = $@"mpsa_{sq.name}_{format}_overall_distribution_{item.group_suffix}",
                                member = $"{sq_index}_{x.ss}_{item.member_suffix}",
                                perspective = "default",
                                feature_value = x.dist_value,
                            };
                            features.Add(dist);

                            var prob_all = new feature_info()
                            {
                                alphabet = "Overall",
                                dimension = 1,
                                category = "mpsa",
                                source = source.ToString(),
                                @group = $@"mpsa_{sq.name}_all_overall_probability_{item.group_suffix}",
                                member = $"{sq_index}_{format}_{x.ss}_{item.member_suffix}",
                                perspective = "default",
                                feature_value = x.prob_value,
                            };
                            features_prob_all.Add(prob_all);

                            var dist_all = new feature_info()
                            {
                                alphabet = "Overall",
                                dimension = 1,
                                category = "mpsa",
                                source = source.ToString(),
                                @group = $@"mpsa_{sq.name}_all_overall_distribution_{item.group_suffix}",
                                member = $"{sq_index}_{format}_{x.ss}_{item.member_suffix}",
                                perspective = "default",
                                feature_value = x.dist_value,
                            };
                            features_dist_all.Add(dist_all);
                        }
                    }

                    var ss_probabilities_per_aa = reader.ss_probabilities_per_aa;
                    //var split_ss_probabilities_per_aa = reader.split_ss_probabilities_per_aa;

                    var merged_ss_probabilities_per_aa = new List<(string group_suffix, string member_suffix, List<(char ss, int alphabet_id, string alphabet_name, string alphabet_group, double prob_value, double dist_value)> data)>();
                    merged_ss_probabilities_per_aa.Add(( /*"unsplit"*/"", "", ss_probabilities_per_aa));
                    //merged_ss_probabilities_per_aa.AddRange(split_ss_probabilities_per_aa.Select((a, i) => ("split", $"{i}", a)).ToList());


                    foreach (var item in merged_ss_probabilities_per_aa)
                    {
                        foreach (var alphabet in ss_probabilities_per_aa.GroupBy(a => a.alphabet_id).ToList())
                        {
                            var alphabet_id = alphabet.Key;

                            var group = alphabet.ToList();

                            var probs = @group.Select(a => new feature_info()
                            {
                                alphabet = a.alphabet_name,
                                dimension = 1,
                                category = "mpsa",
                                source = source.ToString(),
                                @group = $@"mpsa_{sq.name}_{format}_probability_{a.alphabet_name}_{item.group_suffix}",
                                member = $"{sq_index}_{a.ss}_{a.alphabet_group}_{item.member_suffix}",
                                perspective = "default",
                                feature_value = a.prob_value,
                            }).ToList();

                            if (probs.Count <= max_features) features.AddRange(probs);

                            var dists = @group.Select(a => new feature_info()
                            {
                                alphabet = a.alphabet_name,
                                dimension = 1,
                                category = "mpsa",
                                source = source.ToString(),
                                @group = $@"mpsa_{sq.name}_{format}_distribution_{a.alphabet_name}_{item.group_suffix}",
                                member = $"{sq_index}_{a.ss}_{a.alphabet_group}_{item.member_suffix}",
                                perspective = "default",
                                feature_value = a.dist_value,
                            }).ToList();

                            if (dists.Count <= max_features) features.AddRange(dists);


                            // all

                            var probs_all = @group.Select(a => new feature_info()
                            {
                                alphabet = a.alphabet_name,
                                dimension = 1,
                                category = "mpsa",
                                source = source.ToString(),
                                @group = $@"mpsa_{sq.name}_all_probability_{a.alphabet_name}_{item.group_suffix}",
                                member = $"{sq_index}_{format}_{a.ss}_{a.alphabet_group}_{item.member_suffix}",
                                perspective = "default",
                                feature_value = a.prob_value,
                            }).ToList();

                            if (probs_all.Count <= max_features) features_prob_aa_all.AddRange(probs_all);

                            var dists_all = @group.Select(a => new feature_info()
                            {
                                alphabet = a.alphabet_name,
                                dimension = 1,
                                category = "mpsa",
                                source = source.ToString(),
                                @group = $@"mpsa_{sq.name}_all_distribution_{a.alphabet_name}_{item.group_suffix}",
                                member = $"{sq_index}_{format}_{a.ss}_{a.alphabet_group}_{item.member_suffix}",
                                perspective = "default",
                                feature_value = a.dist_value,
                            }).ToList();

                            if (dists_all.Count <= max_features) features_dist_aa_all.AddRange(dists_all);
                        }
                    }
                }


                foreach (var a in features_dist_all.GroupBy(a => (a.@group, a.alphabet)).Select(feature_infos => feature_infos.ToList()).Where(a => a.Count <= max_features))
                {
                    features.AddRange(a);
                }

                foreach (var a in features_prob_all.GroupBy(a => (a.@group, a.alphabet)).Select(feature_infos => feature_infos.ToList()).Where(a => a.Count <= max_features))
                {
                    features.AddRange(a);
                }

                foreach (var a in features_prob_aa_all.GroupBy(a => (a.@group, a.alphabet)).Select(feature_infos => feature_infos.ToList()).Where(a => a.Count <= max_features))
                {
                    features.AddRange(a);
                }

                foreach (var a in features_dist_aa_all.GroupBy(a => (a.@group, a.alphabet)).Select(feature_infos => feature_infos.ToList()).Where(a => a.Count <= max_features))
                {
                    features.AddRange(a);
                }
            }

            if (calculate_mpsa_classification_data_template == null)
            {
                calculate_mpsa_classification_data_template = features.Select(a => new feature_info(a) { source = "", feature_value = 0 }).ToList();
            }

            return features;
        }

        public static List<feature_info> calculate_ring_classification_data_template = null;

        public static List<feature_info> calculate_ring_classification_data(subsequence_classification_data scd, List<Atom> subsequence_master_atoms, protein_data_sources source, int max_features)
        {
#if DEBUG
            if (Program.verbose_debug) Program.WriteLine($"{nameof(calculate_ring_classification_data)}(subsequence_classification_data scd, List<Atom> subsequence_master_atoms, protein_data_sources source, int max_features);");
#endif

            var features = new List<feature_info>();

            if (subsequence_master_atoms == null || subsequence_master_atoms.Count == 0)
            {
                if (calculate_ring_classification_data_template == null) throw new Exception();

                var template = calculate_ring_classification_data_template.Select(a => new feature_info(a)
                {
                    source = source.ToString(),
                    feature_value = 0
                }).ToList();

                return template;
            }

            var edges = subsequence_master_atoms.SelectMany(a => a.monomer_ring_edges).ToList();
            var nodes = subsequence_master_atoms.SelectMany(a => a.monomer_ring_nodes).ToList();


            // node features

            foreach (var alphabet in feature_calcs.aa_alphabets_inc_overall)
            {
                foreach (var alphabet_group in alphabet.groups)
                {
                    var nodes_a = nodes.Where(a => alphabet_group.group_amino_acids.Contains(a.Residue1)).ToList();


                    var feats = new List<feature_info>();


                    var degrees = nodes_a.Select(a => a.Degree).ToArray();
                    var degrees_ds = descriptive_stats.get_stat_values(degrees, nameof(degrees));
                    var degrees_ds_e = descriptive_stats.encode(degrees_ds);
                    var degrees_ds_e_f = degrees_ds_e.Select(a => new feature_info()
                    {
                        alphabet = alphabet.name,
                        dimension = 3,
                        category = $"ring_nodes",
                        source = source.ToString(),
                        group = $"ring_nodes_{nameof(degrees)}_{alphabet.name}",
                        member = a.member_id,
                        perspective = a.perspective_id,
                        feature_value = a.perspective_value
                    }).ToList();
                    feats.AddRange(degrees_ds_e_f);

                    var rapdf = nodes_a.Select(a => a.Rapdf).ToArray();
                    var rapdf_ds = descriptive_stats.get_stat_values(rapdf, nameof(rapdf));
                    var rapdf_ds_e = descriptive_stats.encode(rapdf_ds);
                    var rapdf_ds_e_f = rapdf_ds_e.Select(a => new feature_info()
                    {
                        alphabet = alphabet.name,
                        dimension = 3,
                        category = $"ring_nodes",
                        source = source.ToString(),
                        group = $"ring_nodes_{nameof(rapdf)}_{alphabet.name}",
                        member = a.member_id,
                        perspective = a.perspective_id,
                        feature_value = a.perspective_value
                    }).ToList();
                    feats.AddRange(rapdf_ds_e_f);

                    var tap = nodes_a.Select(a => a.Tap).ToArray();
                    var tap_ds = descriptive_stats.get_stat_values(tap, nameof(tap));
                    var tap_ds_e = descriptive_stats.encode(tap_ds);
                    var tap_ds_e_f = tap_ds_e.Select(a => new feature_info()
                    {
                        alphabet = alphabet.name,
                        dimension = 3,
                        category = $"ring_nodes",
                        source = source.ToString(),
                        group = $"ring_nodes_{nameof(tap)}_{alphabet.name}",
                        member = a.member_id,
                        perspective = a.perspective_id,
                        feature_value = a.perspective_value
                    }).ToList();
                    feats.AddRange(tap_ds_e_f);




                    features.AddRange(feats);

                }
            }



            // edge features

            var dir_types = new[]
            {
                "MC", "SC", "all"
            };
            var bond_types = new[]
            {
                "HBOND",
                "IONIC",
                "PICATION",
                "PIPISTACK",
                "SSBOND",
                "VDW",
                "all"
            };

            foreach (var alphabet in feature_calcs.aa_alphabets_inc_overall)
            {
                foreach (var alphabet_group in alphabet.groups)
                {
                    foreach (var bond_type in bond_types)
                    {
                        var bonds = edges.Where(a => bond_type == "all" || a.Interaction.interaction_type == bond_type)
                            .ToList();

                        bonds = bonds.Where(a => alphabet_group.group_amino_acids.Contains(a.NodeId1.amino_acid1)).ToList();

                        foreach (var dir_type1 in dir_types)
                        {
                            var bonds1 = bonds
                                .Where(a => dir_type1 == "all" || a.Interaction.subtype_node1 == dir_type1)
                                .ToList();

                            foreach (var dir_type2 in dir_types)
                            {
                                var bonds2 = bonds1
                                    .Where(a => dir_type2 == "all" || a.Interaction.subtype_node2 == dir_type2)
                                    .ToList();

                                var count = bonds2.Count;

                                var distances = bonds2.Select(a => a.Distance).ToArray();
                                var distances_ds = descriptive_stats.get_stat_values(distances, nameof(distances));
                                var distnaces_ds_e = descriptive_stats.encode(distances_ds);

                                var angles = bonds2.Where(a => a.Angle != null).Select(a => a.Angle.Value).ToArray();
                                var angles_ds = descriptive_stats.get_stat_values(angles, nameof(angles));
                                var angles_ds_e = descriptive_stats.encode(angles_ds);

                                var energies = bonds2.Select(a => a.Energy).ToArray();
                                var energies_ds = descriptive_stats.get_stat_values(energies, nameof(energies));
                                var energies_ds_e = descriptive_stats.encode(energies_ds);

                                var feats = new List<feature_info>();

                                var count_f = new feature_info()
                                {
                                    alphabet = alphabet.name,
                                    dimension = 3,
                                    category = "ring_edges",
                                    source = source.ToString(),
                                    @group = $"ring_edges_{nameof(count)}_{bond_type}_{dir_type1}_{dir_type2}_{alphabet.name}",
                                    member = "count",
                                    perspective = "default",
                                    feature_value = count
                                };
                                feats.Add(count_f);

                                var distances_ds_e_f = distnaces_ds_e.Select(a => new feature_info()
                                {
                                    alphabet = alphabet.name,
                                    dimension = 3,
                                    category = $"ring_edges",
                                    source = source.ToString(),
                                    group = $"ring_edges_{nameof(distances)}_{bond_type}_{dir_type1}_{dir_type2}_{alphabet.name}",
                                    member = a.member_id,
                                    perspective = a.perspective_id,
                                    feature_value = a.perspective_value
                                }).ToList();
                                feats.AddRange(distances_ds_e_f);


                                var angles_ds_e_f = angles_ds_e.Select(a => new feature_info()
                                {
                                    alphabet = alphabet.name,
                                    dimension = 3,
                                    category = $"ring_edges",
                                    source = source.ToString(),
                                    group = $"ring_edges_{nameof(angles)}_{bond_type}_{dir_type1}_{dir_type2}_{alphabet.name}",
                                    member = a.member_id,
                                    perspective = a.perspective_id,
                                    feature_value = a.perspective_value
                                }).ToList();
                                feats.AddRange(angles_ds_e_f);


                                var energies_ds_e_f = energies_ds_e.Select(a => new feature_info()
                                {
                                    alphabet = alphabet.name,
                                    dimension = 3,
                                    category = $"ring_edges",
                                    source = source.ToString(),
                                    group = $"ring_edges_{nameof(energies)}_{bond_type}_{dir_type1}_{dir_type2}_{alphabet.name}",
                                    member = a.member_id,
                                    perspective = a.perspective_id,
                                    feature_value = a.perspective_value
                                }).ToList();
                                feats.AddRange(energies_ds_e_f);




                                features.AddRange(feats);
                            }
                        }
                    }
                }
            }


            if (calculate_ring_classification_data_template == null)
            {
                calculate_ring_classification_data_template = features.Select(a => new feature_info(a) { source = "", feature_value = 0 }).ToList();
            }

            return features;
        }

        public static List<feature_info> calculate_foldx_classification_data_subsequence_3d_template = null;
        public static List<feature_info> calculate_foldx_classification_data_neighbourhood_3d_template = null;
        public static List<feature_info> calculate_foldx_classification_data_protein_3d_template = null;

        //if (source != protein_data_sources.protein_3d)
        public static List<feature_info> calculate_foldx_classification_data(subsequence_classification_data scd, List<Atom> subsequence_master_atoms, protein_data_sources source, int max_features)
        {
#if DEBUG
            if (Program.verbose_debug) Program.WriteLine($"{nameof(calculate_foldx_classification_data)}(subsequence_classification_data scd, List<Atom> subsequence_master_atoms, protein_data_sources source, int max_features);");
#endif

            if (subsequence_master_atoms == null || subsequence_master_atoms.Count == 0)
            {
                if (source == protein_data_sources.subsequence_3d)
                {
                    if (calculate_foldx_classification_data_subsequence_3d_template == null) throw new Exception();

                    var template = calculate_foldx_classification_data_subsequence_3d_template.Select(a => new feature_info(a)
                    {
                        source = source.ToString(),
                        feature_value = 0
                    }).ToList();

                    return template;
                }
                else if (source == protein_data_sources.neighbourhood_3d)
                {
                    if (calculate_foldx_classification_data_neighbourhood_3d_template == null) throw new Exception();

                    var template = calculate_foldx_classification_data_neighbourhood_3d_template.Select(a => new feature_info(a)
                    {
                        source = source.ToString(),
                        feature_value = 0
                    }).ToList();

                    return template;
                }
                else if (source == protein_data_sources.protein_3d)
                {
                    if (calculate_foldx_classification_data_protein_3d_template == null) throw new Exception();

                    var template = calculate_foldx_classification_data_protein_3d_template.Select(a => new feature_info(a)
                    {
                        source = source.ToString(),
                        feature_value = 0
                    }).ToList();

                    return template;

                }
            }

            var features = new List<feature_info>();

            //var make_protein_foldx_ala_scan_feature = true;
            var make_subsequence_foldx_ala_scan_feature = true;
            var make_subsequence_foldx_position_scan_feature = true;
            var make_subsequence_foldx_buildmodel_position_scan_feature = true;
            var make_subsequence_foldx_buildmodel_subsequence_replacement_feature = true;



            var foldx_energy_differences = scd.foldx_energy_differences;
            var foldx_residues_aa_mutable = foldx_caller.foldx_residues_aa_mutable;


            //var amino_acids = "ARNDCQEGHILKMFPSTWYV";
            //var foldx_amino_acids = string.Join("", foldx_residues_aa_mutable.Select(a => a.foldx_aa_code1).Distinct().ToList());
            //var foldx_specific_amino_acids = string.Join("",foldx_amino_acids.Except(amino_acids).ToList());

            //var alphabets = feature_calcs.aa_alphabets.ToList();
            var aa_alphabets_inc_overall_foldx = feature_calcs.aa_alphabets_inc_overall_foldx.ToList();
            //alphabets.Add((-1, "Overall", new List<string>() { foldx_amino_acids }));
            //alphabets = alphabets.Where(a => !String.Equals(a.name, "Normal", StringComparison.InvariantCultureIgnoreCase)).ToList();
            //alphabets = alphabets.Where(a => a.groups.Count <= 4).ToList();

            //var foldx_alphabets = feature_calcs.aa_alphabets.ToList();
            //foldx_alphabets.Add((-1, "Overall", new List<string>() { foldx_amino_acids }));


            if (make_subsequence_foldx_ala_scan_feature)
            {
                // ALA SCANNING: ALA substitution for each interface amino acid
                var foldx_cmd = $@"foldx_ala_scan";

                var foldx_ala_scan_feats = new List<feature_info>();


                var foldx_ala_scanning_result = foldx_energy_differences.foldx_ala_scanning_result_subsequence.data
                    .OrderBy(a => a.residue_index).ToList();
                var foldx_ala_scanning_result_split = feature_calcs.split_sequence(foldx_ala_scanning_result, 3, 0, false);

                var foldx_ala = foldx_ala_scanning_result_split.Select(a => (name: $@"split", items: a)).ToList();
                foldx_ala.Add((name: $@"unsplit", items: foldx_ala_scanning_result));

                for (var sq_index = 0; sq_index < foldx_ala.Count; sq_index++)
                {
                    var sq = foldx_ala[sq_index];
                    foreach (var alphabet in aa_alphabets_inc_overall_foldx)
                    {
                        foreach (var alphabet_group in alphabet.groups)
                        {
                            // 1. Overall - get ddg of sequence amino acids where the amino acid is ANY amino acid.
                            var items = sq.items.Where(a => alphabet_group.group_amino_acids.Contains(a.original_foldx_amino_acid_1)).ToList();
                            var items_ddg = items.Select(a => a.ddg).ToArray();


                            var items_ddg_ds = descriptive_stats.get_stat_values(items_ddg, $"{alphabet_group.group_name}");
                            var items_ddg_ds_encoded = descriptive_stats.encode(items_ddg_ds);
                            var items_ddg_ds_encoded_features = items_ddg_ds_encoded.Select(a => new feature_info()
                            {
                                alphabet = alphabet.name,
                                dimension = 3,
                                category = $@"{foldx_cmd}",
                                source = source.ToString(),
                                @group = $@"{foldx_cmd}_sequence_{sq.name}_{alphabet.name}",
                                member = $@"{sq_index}_{a.member_id}",
                                perspective = a.perspective_id,
                                feature_value = a.perspective_value
                            }).ToList();
                            foldx_ala_scan_feats.AddRange(items_ddg_ds_encoded_features);
                        }
                    }
                }

                features.AddRange(foldx_ala_scan_feats);

            }



            if (make_subsequence_foldx_position_scan_feature)
            {
                if (source == protein_data_sources.subsequence_3d)
                {
                    // POSITION SCANNING: 20AA+11extras substitution for each interface amino acid (note: ignore positions - focus on amino acids)
                    // feat1: average energy difference overall (All mutation amino acids, any position, whole sequence)
                    // feat2: average energy difference per mutant amino acid type i.e. 31 (All mutation amino acids, specific positions, whole sequence)
                    // feat3: average energy difference per sequence amino acid i.e. 20 

                    var foldx_cmd = $@"foldx_pos_scan";

                    var foldx_pos_feats = new List<feature_info>();


                    var foldx_pos_scanning_result = foldx_energy_differences.foldx_position_scanning_result_subsequence.data.OrderBy(a => a.residue_index).ToList();
                    var foldx_pos_scanning_result_split = feature_calcs.split_sequence(foldx_pos_scanning_result, 3, 0, false);

                    var foldx_pos = foldx_pos_scanning_result_split.Select(a => (name: $@"split", items: a)).ToList();
                    foldx_pos.Add((name: $@"unsplit", items: foldx_pos_scanning_result));

                    for (var sq_index = 0; sq_index < foldx_pos.Count; sq_index++)
                    {
                        var sq = foldx_pos[sq_index];
                        foreach (var alphabet in aa_alphabets_inc_overall_foldx)
                        {
                            // todo: consider looping alphabet for both row/column, problem: it makes a lot of extra features, however: the features would be more specific i.e. sensitive

                            foreach (var g_row_original_foldx_amino_acid_alphabet_group in alphabet.groups)
                            {
                                // rows: compare by original amino acids (compact seq of len L to 20..31 AA) note: position scan doesn't stick to the 20 standard AA
                                //if (foldx_specific_amino_acids.All(a => "" + a != g_row_original_amino_acid))
                                {
                                    // 1. Overall - get ddg of sequence amino acids where the amino acid is ANY amino acid.
                                    var items = sq.items.Where(a => g_row_original_foldx_amino_acid_alphabet_group.group_amino_acids.Contains(a.original_foldx_amino_acid_1)).ToList();
                                    var items_ddg = items.Select(a => a.ddg).ToArray();


                                    var items_ddg_ds = descriptive_stats.get_stat_values(items_ddg, $"{g_row_original_foldx_amino_acid_alphabet_group.group_name}");
                                    var items_ddg_ds_encoded = descriptive_stats.encode(items_ddg_ds);
                                    var items_ddg_ds_encoded_features = items_ddg_ds_encoded.Select(a => new feature_info()
                                    {
                                        alphabet = alphabet.name,
                                        dimension = 3,
                                        category = $@"{foldx_cmd}",
                                        source = source.ToString(),
                                        @group = $@"{foldx_cmd}_sequence_{sq.name}_{alphabet.name}",
                                        member = $@"{sq_index}_{a.member_id}",
                                        perspective = a.perspective_id,
                                        feature_value = a.perspective_value
                                    }).ToList();
                                    foldx_pos_feats.AddRange(items_ddg_ds_encoded_features);
                                }

                                // columns: compare by mutation amino acids
                                {
                                    var g_col_foldx_amino_acid_alphabet_group = g_row_original_foldx_amino_acid_alphabet_group;

                                    var items = sq.items.Where(a => g_col_foldx_amino_acid_alphabet_group.group_amino_acids.Contains(a.mutant_foldx_amino_acid_1)).ToList();
                                    var items_ddg = items.Select(a => a.ddg).ToArray();


                                    var items_ddg_ds = descriptive_stats.get_stat_values(items_ddg, $"{g_col_foldx_amino_acid_alphabet_group.group_name}");
                                    var items_ddg_ds_encoded = descriptive_stats.encode(items_ddg_ds);
                                    var items_ddg_ds_encoded_features = items_ddg_ds_encoded.Select(a => new feature_info()
                                    {
                                        alphabet = alphabet.name,
                                        dimension = 3,
                                        category = $@"{foldx_cmd}",
                                        source = source.ToString(),
                                        @group = $@"{foldx_cmd}_mutants_{sq.name}_{alphabet.name}",
                                        member = $@"{sq_index}_{a.member_id}",
                                        perspective = a.perspective_id,
                                        feature_value = a.perspective_value
                                    }).ToList();
                                    foldx_pos_feats.AddRange(items_ddg_ds_encoded_features);
                                }

                                // matrix row:column code: compare by cross-reference of original amino acids and mutation amino acids

                                foreach (var g_col_mutant_foldx_amino_acid_alphabet_group in alphabet.groups)
                                {
                                    var items = sq.items.Where(a => g_row_original_foldx_amino_acid_alphabet_group.group_amino_acids.Contains(a.original_foldx_amino_acid_1) && g_col_mutant_foldx_amino_acid_alphabet_group.group_amino_acids.Contains(a.mutant_foldx_amino_acid_1)).ToList();
                                    var items_ddg = items.Select(a => a.ddg).ToArray();


                                    var items_ddg_ds = descriptive_stats.get_stat_values(items_ddg, $"{g_row_original_foldx_amino_acid_alphabet_group}_{g_col_mutant_foldx_amino_acid_alphabet_group.group_name}");
                                    var items_ddg_ds_encoded = descriptive_stats.encode(items_ddg_ds);
                                    var items_ddg_ds_encoded_features = items_ddg_ds_encoded.Select(a => new feature_info()
                                    {
                                        alphabet = alphabet.name,
                                        dimension = 3,
                                        category = $"{foldx_cmd}",
                                        source = source.ToString(),
                                        @group = $"{foldx_cmd}_{sq.name}_{alphabet.name}",
                                        member = $@"{sq_index}_{a.member_id}",
                                        perspective = a.perspective_id,
                                        feature_value = a.perspective_value
                                    }).ToList();
                                    foldx_pos_feats.AddRange(items_ddg_ds_encoded_features);
                                }
                            }
                        }

                        features.AddRange(foldx_pos_feats);
                    }
                }
            }

            if (make_subsequence_foldx_buildmodel_position_scan_feature)
            {
                if (source == protein_data_sources.subsequence_3d)
                {

                    // divide by 31 (number of foldx amino acid codes), spli

                    var foldx_cmd = $@"foldx_bm_ps";
                    var foldx_bm_ps_feats = new List<feature_info>();


                    var foldx_bm_ps_scanning_result = foldx_energy_differences.foldx_buildmodel_position_scan_result_subsequence.data.OrderBy(a => a.mutation_positions_data.residue_index).ToList();


                    var foldx_bm_ps_scanning_result_split_grouped = feature_calcs.split_sequence(foldx_bm_ps_scanning_result.GroupBy(a => a.mutation_positions_data.residue_index).Select(a => a.ToList()).ToList(), 3, 0, false);

                    var foldx_bm_ps_scanning_result_split = foldx_bm_ps_scanning_result_split_grouped.Select(a => a.SelectMany(b => b).ToList()).ToList();

                    var foldx_bm_ps = foldx_bm_ps_scanning_result_split.Select(a => (name: $@"split", items: a)).ToList();

                    foldx_bm_ps.Add((name: $@"unsplit", items: foldx_bm_ps_scanning_result));


                    for (var sq_index = 0; sq_index < foldx_bm_ps.Count; sq_index++)
                    {
                        var sq = foldx_bm_ps[sq_index];
                        foreach (var alphabet in aa_alphabets_inc_overall_foldx)
                        {
                            // todo: consider looping alphabet for both row/column, problem: it makes a lot of extra features, however: the features would be more specific i.e. sensitive

                            foreach (var g_row_original_foldx_amino_acid_alphabet_group in alphabet.groups)
                            {
                                // rows: compare by original amino acids (compact seq of len L to 20..31 AA) note: position scan doesn't stick to the 20 standard AA
                                //if (foldx_specific_amino_acids.All(a => "" + a != g_row_original_amino_acid))
                                {
                                    // 1. Overall - get ddg of sequence amino acids where the amino acid is ANY amino acid.
                                    var items = sq.items.Where(a => g_row_original_foldx_amino_acid_alphabet_group.group_amino_acids.Contains(a.mutation_positions_data.original_amino_acid1)).ToList();
                                    var items_ddg_list = items.SelectMany(a => a.properties()).GroupBy(a => a.name).Select(a => (energy_name: a.Key, values: a.Select(b => b.value).ToArray())).ToList();

                                    foreach (var items_ddg in items_ddg_list)
                                    {
                                        var items_ddg_ds = descriptive_stats.get_stat_values(items_ddg.values, $"{g_row_original_foldx_amino_acid_alphabet_group.group_name}");
                                        var items_ddg_ds_encoded = descriptive_stats.encode(items_ddg_ds);
                                        var items_ddg_ds_encoded_features_separated = items_ddg_ds_encoded.Select(a => new feature_info()
                                        {
                                            alphabet = alphabet.name,
                                            dimension = 3,
                                            category = $@"{foldx_cmd}",
                                            source = source.ToString(),
                                            @group = $@"{foldx_cmd}_sequence_{items_ddg.energy_name}_{sq.name}_{alphabet.name}",
                                            member = $@"{sq_index}_{items_ddg.energy_name}_{a.member_id}",
                                            perspective = a.perspective_id,
                                            feature_value = a.perspective_value
                                        }).ToList();
                                        foldx_bm_ps_feats.AddRange(items_ddg_ds_encoded_features_separated);

                                        var items_ddg_ds_encoded_features = items_ddg_ds_encoded.Select(a => new feature_info()
                                        {
                                            alphabet = alphabet.name,
                                            dimension = 3,
                                            category = $@"{foldx_cmd}",
                                            source = source.ToString(),
                                            @group = $@"{foldx_cmd}_sequence_{sq.name}_{alphabet.name}",
                                            member = $@"{sq_index}_{items_ddg.energy_name}_{a.member_id}",
                                            perspective = a.perspective_id,
                                            feature_value = a.perspective_value
                                        }).ToList();
                                        foldx_bm_ps_feats.AddRange(items_ddg_ds_encoded_features);
                                    }
                                }

                                // columns: compare by mutation amino acids
                                {
                                    var g_col_foldx_amino_acid = g_row_original_foldx_amino_acid_alphabet_group;

                                    var items = sq.items.Where(a => g_col_foldx_amino_acid.group_amino_acids.Contains(a.mutation_positions_data.mutant_foldx_amino_acid1)).ToList();
                                    var items_ddg_list = items.SelectMany(a => a.properties()).GroupBy(a => a.name).Select(a => (energy_name: a.Key, values: a.Select(b => b.value).ToArray())).ToList();

                                    foreach (var items_ddg in items_ddg_list)
                                    {
                                        var items_ddg_ds = descriptive_stats.get_stat_values(items_ddg.values, $"{g_col_foldx_amino_acid}");
                                        var items_ddg_ds_encoded = descriptive_stats.encode(items_ddg_ds);
                                        var items_ddg_ds_encoded_features_separated = items_ddg_ds_encoded.Select(a => new feature_info()
                                        {
                                            alphabet = alphabet.name,
                                            dimension = 3,
                                            category = $@"{foldx_cmd}",
                                            source = source.ToString(),
                                            @group = $@"{foldx_cmd}_mutants_{items_ddg.energy_name}_{sq.name}_{alphabet.name}",
                                            member = $@"{sq_index}_{items_ddg.energy_name}_{a.member_id}",
                                            perspective = a.perspective_id,
                                            feature_value = a.perspective_value
                                        }).ToList();
                                        foldx_bm_ps_feats.AddRange(items_ddg_ds_encoded_features_separated);

                                        var items_ddg_ds_encoded_features = items_ddg_ds_encoded.Select(a => new feature_info()
                                        {
                                            alphabet = alphabet.name,
                                            dimension = 3,
                                            category = $@"{foldx_cmd}",
                                            source = source.ToString(),
                                            @group = $@"{foldx_cmd}_mutants_{sq.name}_{alphabet.name}",
                                            member = $@"{sq_index}_{items_ddg.energy_name}_{a.member_id}",
                                            perspective = a.perspective_id,
                                            feature_value = a.perspective_value
                                        }).ToList();
                                        foldx_bm_ps_feats.AddRange(items_ddg_ds_encoded_features);
                                    }
                                }

                                // matrix row:column code: compare by cross-reference of original amino acids and mutation amino acids

                                foreach (var g_col_mutant_foldx_amino_acid_alphabet_group in alphabet.groups)
                                {
                                    var items = sq.items.Where(a => g_row_original_foldx_amino_acid_alphabet_group.group_amino_acids.Contains(a.mutation_positions_data.original_amino_acid1) && g_col_mutant_foldx_amino_acid_alphabet_group.group_amino_acids.Contains(a.mutation_positions_data.mutant_foldx_amino_acid1)).ToList();
                                    var items_ddg_list = items.SelectMany(a => a.properties()).GroupBy(a => a.name).Select(a => (energy_name: a.Key, values: a.Select(b => b.value).ToArray())).ToList();


                                    foreach (var items_ddg in items_ddg_list)
                                    {
                                        var items_ddg_ds = descriptive_stats.get_stat_values(items_ddg.values, $"{g_row_original_foldx_amino_acid_alphabet_group}_{g_col_mutant_foldx_amino_acid_alphabet_group.group_name}");
                                        var items_ddg_ds_encoded = descriptive_stats.encode(items_ddg_ds);
                                        var items_ddg_ds_encoded_features_separated = items_ddg_ds_encoded.Select(a => new feature_info()
                                        {
                                            alphabet = alphabet.name,
                                            dimension = 3,
                                            category = $@"{foldx_cmd}",
                                            source = source.ToString(),
                                            @group = $@"{foldx_cmd}_{items_ddg.energy_name}_{sq.name}_{alphabet.name}",
                                            member = $@"{sq_index}_{items_ddg.energy_name}_{a.member_id}",
                                            perspective = a.perspective_id,
                                            feature_value = a.perspective_value
                                        }).ToList();
                                        foldx_bm_ps_feats.AddRange(items_ddg_ds_encoded_features_separated);

                                        var items_ddg_ds_encoded_features = items_ddg_ds_encoded.Select(a => new feature_info()
                                        {
                                            alphabet = alphabet.name,
                                            dimension = 3,
                                            category = $@"{foldx_cmd}",
                                            source = source.ToString(),
                                            @group = $@"{foldx_cmd}_{sq.name}_{alphabet.name}",
                                            member = $@"{sq_index}_{items_ddg.energy_name}_{a.member_id}",
                                            perspective = a.perspective_id,
                                            feature_value = a.perspective_value
                                        }).ToList();

                                        foldx_bm_ps_feats.AddRange(items_ddg_ds_encoded_features);
                                    }
                                }
                            }
                        }

                        features.AddRange(foldx_bm_ps_feats);
                    }
                }
            }

            if (make_subsequence_foldx_buildmodel_subsequence_replacement_feature)
            {

                if (source == protein_data_sources.subsequence_3d)
                {
                    // divide by 31 (number of foldx amino acid codes), split

                    // note: there is no way to split this data - it is 1 data point per interface

                    var foldx_cmd = "foldx_bm_if_sub";
                    var foldx_bm_if_sub_feats = new List<feature_info>();


                    var foldx_bm_if_sub = foldx_energy_differences.foldx_buildmodel_subsequence_mutant_result_subsequence.data.ToList();





                    foreach (var alphabet in aa_alphabets_inc_overall_foldx)
                    {
                        // todo: consider looping alphabet for both row/column, problem: it makes a lot of extra features, however: the features would be more specific i.e. sensitive

                        foreach (var g_row_mutant_foldx_amino_acid_alphabet_group in alphabet.groups)
                        {
                            // rows: compare by original amino acids (compact seq of len L to 20..31 AA) note: position scan doesn't stick to the 20 standard AA
                            //if (foldx_specific_amino_acids.All(a => "" + a != g_row_original_amino_acid))
                            {
                                // 1. Overall - get ddg of sequence amino acids where the amino acid is ANY amino acid.
                                var items = foldx_bm_if_sub.Where(a => g_row_mutant_foldx_amino_acid_alphabet_group.group_amino_acids.Contains(a.mutation_positions_data.First().mutant_foldx_amino_acid1)).ToList();
                                var items_ddg_list = items.SelectMany(a => a.properties()).GroupBy(a => a.name).Select(a => (energy_name: a.Key, values: a.Select(b => b.value).ToArray())).ToList();

                                foreach (var items_ddg in items_ddg_list)
                                {


                                    var items_ddg_ds = descriptive_stats.get_stat_values(items_ddg.values, $"{g_row_mutant_foldx_amino_acid_alphabet_group.group_name}");
                                    var items_ddg_ds_encoded = descriptive_stats.encode(items_ddg_ds);
                                    var items_ddg_ds_encoded_features_separated = items_ddg_ds_encoded.Select(a => new feature_info()
                                    {
                                        alphabet = alphabet.name,
                                        dimension = 3,
                                        category = $"{foldx_cmd}",
                                        source = source.ToString(),
                                        group = $"{foldx_cmd}_mutant_{items_ddg.energy_name}_{alphabet.name}",
                                        member = $"{items_ddg.energy_name}_{a.member_id}",
                                        perspective = a.perspective_id,
                                        feature_value = a.perspective_value
                                    }).ToList();
                                    foldx_bm_if_sub_feats.AddRange(items_ddg_ds_encoded_features_separated);

                                    var items_ddg_ds_encoded_features = items_ddg_ds_encoded.Select(a => new feature_info()
                                    {
                                        alphabet = alphabet.name,
                                        dimension = 3,
                                        category = $"{foldx_cmd}",
                                        source = source.ToString(),
                                        group = $"{foldx_cmd}_mutant_{alphabet.name}",
                                        member = $"{items_ddg.energy_name}_{a.member_id}",
                                        perspective = a.perspective_id,
                                        feature_value = a.perspective_value
                                    }).ToList();
                                    foldx_bm_if_sub_feats.AddRange(items_ddg_ds_encoded_features);
                                }
                            }
                        }
                    }

                    features.AddRange(foldx_bm_if_sub_feats);
                }
            }


            if (source == protein_data_sources.subsequence_3d)
            {
                if (calculate_foldx_classification_data_subsequence_3d_template == null)
                {
                    calculate_foldx_classification_data_subsequence_3d_template = features.Select(a => new feature_info(a) { source = "", feature_value = 0 }).ToList();
                }
            }
            else if (source == protein_data_sources.neighbourhood_3d)
            {
                if (calculate_foldx_classification_data_neighbourhood_3d_template == null)
                {
                    calculate_foldx_classification_data_neighbourhood_3d_template = features.Select(a => new feature_info(a) { source = "", feature_value = 0 }).ToList();
                }
            }
            else if (source == protein_data_sources.protein_3d)
            {
                if (calculate_foldx_classification_data_protein_3d_template == null)
                {
                    calculate_foldx_classification_data_protein_3d_template = features.Select(a => new feature_info(a) { source = "", feature_value = 0 }).ToList();
                }
            }

            return features;
        }

        public static List<feature_info> calculate_sequence_geometry_classification_data_template = null;
        public static List<feature_info> calculate_sequence_geometry_classification_data(subsequence_classification_data scd, List<Atom> subsequence_master_atoms, protein_data_sources source, int max_features)
        {
#if DEBUG
            if (Program.verbose_debug) Program.WriteLine($"{nameof(calculate_sequence_geometry_classification_data)}(subsequence_classification_data scd, List<Atom> subsequence_master_atoms, protein_data_sources source, int max_features);");
#endif

            // note: cannot use PDB length data, because that leaks information.

            if (subsequence_master_atoms == null || subsequence_master_atoms.Count == 0)
            {
                if (calculate_sequence_geometry_classification_data_template == null) throw new Exception();

                var template = calculate_sequence_geometry_classification_data_template.Select(a => new feature_info(a)
                {
                    source = source.ToString(),
                    feature_value = 0
                }).ToList();

                return template;
            }

            var features = new List<feature_info>();

            var make_seq_len_feature = true;
            var use_pdb_seq = true;
            var use_uniprot_seq = false;

            //var make_parity_feature = true;
            //var make_subsequence_sequence_position_feature = true;
            //var make_protein_len_feature = true;


            var subsequence = scd?.aa_subsequence ?? "";
            //var uniprot_sequence = use_uniprot_seq ? (scd?.pdb_chain_master_atoms?.FirstOrDefault()?.uniprot_sequence ?? "") : "";
            var pdb_sequence = string.Join("", scd?.pdb_chain_master_atoms?.Select(a => a.amino_acid).ToList() ?? new List<char>());

            //var pdb_id_simple = scd.pdb_id.Substring(0, 4) + scd.chain_id;

            //var uniprot_mapping = File.ReadAllLines($@"C:\betastrands_dataset\uniprot_pdb_mapping\{pdb_id_simple}.txt").Select(a =>
            //{
            //    var s = a.Split(new char[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);
            //    var i = 0;

            //    var x = (
            //            PDB_code: s[i++],
            //            PDB_Chain_label: s[i++][0],
            //            PDB_Sequential_residue_number: int.Parse(s[i++]),
            //            PDB_Residue_type: s[i++],
            //            PDB_Residue_number: int.Parse(s[i++]),
            //            UniProt_Accession: s[i++],
            //            UniProt_Residue_type: s[i++][0],
            //            UniProt_Sequential_residue_number: int.Parse(s[i++])
            //        );

            //    return x;
            //}).ToList();

            var subseq_len = (double)subsequence.Length;
            var pdb_len = (double)pdb_sequence.Length;
            //var uniprot_len = (double)uniprot_sequence.Length;

            var subseq_len_relative_to_pdb_len = (double)(pdb_len == 0 ? 0d : (double)subseq_len / (double)pdb_len);
            //var subseq_len_relative_to_uniprot = (double)(uniprot_len == 0 ? 0d : (double)subseq_len / (double)uniprot_len);

            //var protein_len_rel = protein_len / subseq_len_rel;
            var middle_subseq_res = scd.res_ids[scd.res_ids.Count / 2];

            //var middle_uniprot_index = uniprot_mapping.FirstOrDefault(a => a.PDB_Residue_number == middle_subseq_res.residue_index).UniProt_Sequential_residue_number;
            //var middle_subseq_res_uniprot_index_pct = (middle_uniprot_index == -1 || scd.pdb_chain_master_atoms == null) ? 0d : (double)middle_uniprot_index / (double)uniprot_len;

            var middle_subseq_res_pdb_index = (scd.pdb_chain_master_atoms == null || scd.pdb_chain_master_atoms.Count == 0) ? 0d : scd.pdb_chain_master_atoms.FindIndex(a => a.residue_index == middle_subseq_res.residue_index);

            var middle_subseq_res_pdb_index_pct = (middle_subseq_res_pdb_index == -1 || scd.pdb_chain_master_atoms == null) ? 0d : (double)middle_subseq_res_pdb_index / (double)pdb_len;





            //var length_rel = (string.IsNullOrWhiteSpace(scd.aa_subsequence) ||
            //                  scd.pdb_chain_master_atoms == null ||
            //                  scd.pdb_chain_master_atoms.Count == 0)
            //    ?
            //    0d
            //    :
            //    (double)subseq_len / (double)protein_len;

           


            if (make_seq_len_feature)
            {
                var length_features = new List<feature_info>();

                var category = "length_sequence";
                var alphabet = "Overall";

                var x1 = new feature_info()
                {
                    alphabet = alphabet,
                    dimension = 1,
                    category = category,
                    source = source.ToString(),
                    group = "length_subsequence_abs",
                    member = "length_subsequence_abs",
                    perspective = "default",
                    feature_value = subseq_len
                };

                length_features.Add(x1);


                if (use_pdb_seq)
                {
                    var x2a = new feature_info()
                    {
                        alphabet = alphabet,
                        dimension = 1,
                        category = category,
                        source = source.ToString(),
                        group = "length_subsequence_rel_pdb",
                        member = "length_subsequence_rel_pdb",
                        perspective = "default",
                        feature_value = subseq_len_relative_to_pdb_len
                    };
                    length_features.Add(x2a);
                }

                if (use_uniprot_seq)
                {
                    //var x2b = new feature_info()
                    //{
                    //    alphabet = alphabet,
                    //    dimension = 1,
                    //    category = category,
                    //    source = source.ToString(),
                    //    group = "length_subsequence_rel_uniprot",
                    //    member = "length_subsequence_rel_uniprot",
                    //    perspective = "default",
                    //    feature_value = subseq_len_relative_to_uniprot
                    //};
                    //length_features.Add(x2b);
                }

                if (use_pdb_seq)
                {
                    var x3a = new feature_info()
                    {
                        alphabet = alphabet,
                        dimension = 1,
                        category = category,
                        source = source.ToString(),
                        group = "length_pdb",
                        member = "length_pdb",
                        perspective = "default",
                        feature_value = pdb_len
                    };
                    length_features.Add(x3a);
                }


                if (use_uniprot_seq)
                {
                    //var x3b = new feature_info()
                    //{
                    //    alphabet = alphabet,
                    //    dimension = 1,
                    //    category = category,
                    //    source = source.ToString(),
                    //    group = "length_uniprot",
                    //    member = "length_uniprot",
                    //    perspective = "default",
                    //    feature_value = uniprot_len
                    //};
                    //length_features.Add(x3b);

                }

                if (use_pdb_seq)
                {
                    var x4a = new feature_info()
                    {
                        alphabet = alphabet,
                        dimension = 1,
                        category = category,
                        source = source.ToString(),
                        group = "sequence_position_rel_pdb",
                        member = "sequence_position_rel_pdb",
                        perspective = "default",
                        feature_value = middle_subseq_res_pdb_index_pct
                    };
                    length_features.Add(x4a);
                }


                if (use_uniprot_seq)
                {
                    //var x4b = new feature_info()
                    //{
                    //    alphabet = alphabet,
                    //    dimension = 1,
                    //    category = category,
                    //    source = source.ToString(),
                    //    group = "sequence_position_rel_uniprot",
                    //    member = "default",
                    //    perspective = "default",
                    //    feature_value = middle_subseq_res_uniprot_index_pct
                    //};

                    //length_features.Add(x4b);
                }









                var all_length_features = length_features.Select(a => new feature_info(a)
                {
                    @group = "length_all"
                }).ToList();

                features.AddRange(length_features);
                features.AddRange(all_length_features);


            }

            if (calculate_sequence_geometry_classification_data_template == null)
            {
                calculate_sequence_geometry_classification_data_template = features.Select(a => new feature_info(a) { source = "", feature_value = 0 }).ToList();
            }


            return features;
        }

        public static List<(string name, List<aaindex.aaindex_entry> list)> aaindex_subset_templates = aaindex_subset_templates_search();

        public static List<(string name, List<aaindex.aaindex_entry> list)> aaindex_subset_templates_search()
        {
            var keywords = new List<(string name, string[] list)>();

            // energy terms
            keywords.Add(("energy", new[] { "thermodynamic", "thermodynamics", "thermal", "energy", "gibbs", "solvation", "entropy", "entropies", "energies", "pka", "pk", "ph", "heat", "temperature", "dg", "ddg", "delta-g", "delta g" }));
            // charge terms
            keywords.Add(("charge", new[] { "charge", "polarity", "polar", "charged", "positive", "negative", "electric", "electricity", "electrostatic" }));
            // interaction terms
            keywords.Add(("interaction", new[] { "interaction", "interactions", "attraction", "affinity", "contact", "contacts", "complex", "complexation", "bind", "bond", "bonding", "binding", "bonded", "partner", "partnered", "partnering", "interaction", "intramolecular", "intermolecular", "vdw", "van der waals", "electrostatic", "statics", "hydrogen", "hbond" }));
            // burial/exposed keywords
            keywords.Add(("accessibility", new[] { "buried", "burial", "exposed", "exposure", "hidden", "accessibility", "accessible", "surface", "surfacial", "solvation", "solvent" }));
            // unordered regions keywords
            keywords.Add(("disorder", new[] { "unordered", "disorder", "randomness", "random coil", "random region", "coil", "terminal", "ambiguous", "conformational change" }));
            // beta-strand keywords
            keywords.Add(("strand", new[] { "strand", "sheet", "beta-strand", "beta-sheet", "strand-strand", "sheet-sheet" }));
            // coil keywords
            keywords.Add(("coil", new[] { "coil", "random coil", "unstructured", "unordered", "coiled coil", "terminal coil", "coil-coil", "coil-strand", "coil-helix", "starnd-coil", "helix-coil" }));
            // all secondary structure keywords
            keywords.Add(("ss", new[] { "transformation", "conversion", "conformation", "structure", "helix", "helice", "helical", "coil", "coiled", "helix", "strand", "sheet", "ss", "sec struct", "secondary structure" }));
            // hydrophobicity keywords
            keywords.Add(("hydrophocity", new[] { "hydropathy", "hydrophobe", "hydrophilathy", "hydrophobicity", "hydrophobic", "hydrophil", "hydrophile", "hydrophilic", "hydrophicility" }));
            // composition keywords
            keywords.Add(("composition", new[] { "composition", "propensity", "distribution", "frequency" }));

            List<(string name, List<aaindex.aaindex_entry> list)> aaindices_subsections = keywords.Select(a => 
                (name: a.name,
                    list: 
                    aaindex.aaindex_entries.Where(b => a.list.Any(c => b.D_Data_Description.ToLowerInvariant().Contains(c.ToLowerInvariant()) || b.T_Title_Of_Article.ToLowerInvariant().Contains(c.ToLowerInvariant()))).Distinct().ToList())).Distinct().ToList();

            
            aaindices_subsections.Add(("all", aaindex.aaindex_entries));
            

            // from a paper... which one?
            var dna_binding = new string[] {
            "CHOP780202", "GEIM800106", "PALJ810107", "ZIMJ680104",
            "CIDH920103", "KANM800102", "QIAN880123", "AURR980120",
            "CIDH920105", "KLEP840101", "RACS770103", "MUNV940103",
            "FAUJ880109", "KRIW710101", "RADA880108", "NADH010104",
            "FAUJ880111", "LIFS790101", "ROSM880102", "NADH010106",
            "FINA910104", "MEEJ800101", "SWER830101", "GUYH850105",
            "GEIM800104", "OOBM770102", "ZIMJ680102", "MIYS990104"};
            aaindices_subsections.Add((nameof(dna_binding), aaindex.aaindex_entries.Where(a => dna_binding.Contains(a.H_Accession_Number)).ToList()));

            // from a paper... which one?
            var zernike = new string[] { "BLAM930101", "BIOV880101", "MAXF760101", "TSAJ990101", "NAKH920108", "CEDJ970104", "LIFS790101", "MIYS990104", };
            aaindices_subsections.Add((nameof(zernike), aaindex.aaindex_entries.Where(a => zernike.Contains(a.H_Accession_Number)).ToList()));

            //An Ensemble Method for Predicting Subnuclear Localizations from Primary Protein Structures
            var subnuclear = new string[] { "BULH740101", "BULH740102", "PONP800106", "PONP800104", "PONP800105", "PONP800106", "MANP780101", "EISD840101", "JOND750101", "HOPT810101", "PARJ860101", "JANJ780101", "PONP800107", "CHOC760102", "ROSG850101", "ROSG850102", "BHAR880101", "KARP850101", "KARP850102", "KARP850103", "JANJ780102", "JANJ780103", "LEVM780101", "LEVM780102", "LEVM780103", "GRAR740102", "GRAR740103", "MCMT640101", "PONP800108", "KYTJ820101", };
            aaindices_subsections.Add((nameof(subnuclear), aaindex.aaindex_entries.Where(a => subnuclear.Contains(a.H_Accession_Number)).ToList()));

            //Identification of properties important to protein aggregation using feature selection
            var aggregation = new string[] { "CASG920101", "GUYH850101", "LEVM780102", "PALJ810111", "PONP800105", "PONP800107", "PRAM820103", "PRAM900103", "RICJ880117", "ROBB760110", "ROSM880105", "ROSM880105", "VENT840101", "VHEG790101", "WILM950102", "ZIMJ680101", };
            aaindices_subsections.Add((nameof(aggregation), aaindex.aaindex_entries.Where(a => aggregation.Contains(a.H_Accession_Number)).ToList()));

            //Prediction of Protein–Protein Interaction with Pairwise Kernel Support Vector Machine
            var ppi = new string[] { "LEWP710101", "QIAN880138", "NADH010104", "NAGK730103", "AURR980116" };
            aaindices_subsections.Add((nameof(ppi), aaindex.aaindex_entries.Where(a => ppi.Contains(a.H_Accession_Number)).ToList()));

            //Characterizing informative sequence descriptors and predicting binding affinities of heterodimeric protein complexes
            var affinity = new string[] { "GUYH850105", "SNEP660104", "RACS820113", "MITS020101", "MAXF760103", "CIDH920104", "AURR980119", "TANS770103", "CHOP780101", "PALJ810107", "QIAN880116", "PALJ810110", "TAKK010101", };
            aaindices_subsections.Add((nameof(affinity), aaindex.aaindex_entries.Where(a => affinity.Contains(a.H_Accession_Number)).ToList()));

            // Intersection of the values from papers
            var intersection = new string[] { "PALJ810107", "NADH010104", "LIFS790101", "GUYH850105", "MIYS990104", "PONP800106", "PONP800105", "PONP800107", "LEVM780102", "ROSM880105", };
            aaindices_subsections.Add((nameof(intersection), aaindex.aaindex_entries.Where(a => intersection.Contains(a.H_Accession_Number)).ToList()));


            return aaindices_subsections;
        }

        public static List<feature_info> calculate_aa_index_classification_data_template = null;

        public static List<feature_info> calculate_aa_index_classification_data(List<Atom> subsequence_master_atoms, protein_data_sources source, int max_features)
        {
#if DEBUG
            if (Program.verbose_debug) Program.WriteLine($"{nameof(calculate_aa_index_classification_data)}(List<Atom> subsequence_master_atoms, protein_data_sources source, int max_features);");
#endif


            if (subsequence_master_atoms == null || subsequence_master_atoms.Count == 0)
            {
                if (calculate_aa_index_classification_data_template == null) throw new Exception();

                var template = calculate_aa_index_classification_data_template.Select(a => new feature_info(a)
                {
                    source = source.ToString(),
                    feature_value = 0
                }).ToList();

                return template;
            }

            var aaindices_subsections = aaindex_subset_templates;

            var complete_sequence = string.Join("", subsequence_master_atoms.Select(a => a.amino_acid).ToList());
            var split_sequence = feature_calcs.split_sequence(complete_sequence, 3, false);
            var sequences = new List<(string name, string sequence)>();
            sequences.Add(($@"unsplit", complete_sequence));
            sequences.AddRange(split_sequence.Select(a => ($@"split", a)).ToList());

            var alphabets = feature_calcs.aa_alphabets_inc_overall.ToList();
            //alphabets.Add((-1, "Overall", new List<string>() { $@"ARNDCQEGHILKMFPSTWYV" }));
            alphabets = alphabets.Where(a => !String.Equals(a.name, $@"Normal", StringComparison.InvariantCultureIgnoreCase)).ToList();
            alphabets = alphabets.Where(a => a.groups.Count <= 4).ToList();


            var result = new List<feature_info>();

            // loop through each subsection of the aaindex.  note: one of the subsections is the full aaindex.
            foreach (var aaindexes in aaindices_subsections)
            {
                // loop through the split and unsplit sequences
                for (var sq_index = 0; sq_index < sequences.Count; sq_index++)
                {
                    var sq = sequences[sq_index];
// loop through each aaindex entry without the aaindex subsection
                    foreach (var aaindex_entry in aaindexes.list) //aaindex.aaindex_entries)
                    {
                        // get the values for the sequence 
                        var seq_aaindex_values = aaindex.sequence_aaindex_entry(aaindex_entry.H_Accession_Number, sq.sequence);

                        foreach (var alphabet in alphabets)
                        {
                            foreach (var alphabet_group in alphabet.groups)
                            {
                                var seq_aaindex_values_limited = seq_aaindex_values.Where(a => alphabet_group.group_amino_acids.Contains(a.amino_acid)).Select(a => a.value).ToArray();


                                var ds_values = descriptive_stats.get_stat_values(seq_aaindex_values_limited, "");

                                var e_ds_values = descriptive_stats.encode(ds_values);


                                if (string.Equals(aaindexes.name, "all", StringComparison.InvariantCultureIgnoreCase))
                                {
                                    var f_e_ds_values1 = e_ds_values.Select(a => new feature_info()
                                    {
                                        alphabet = alphabet.name,
                                        dimension = 1,
                                        category = $@"aaindex_{aaindexes.name}",
                                        source = source.ToString(),
                                        @group = $@"aaindex_{aaindexes.name}_{sq.name}_{aaindex_entry.H_Accession_Number}_{alphabet.name}",
                                        member = $@"{sq_index}_{alphabet_group.group_name}_{a.member_id}_{aaindex_entry.H_Accession_Number}",
                                        perspective = a.perspective_id,
                                        feature_value = a.perspective_value
                                    }).ToList();

                                    if (f_e_ds_values1.Count <= max_features) result.AddRange(f_e_ds_values1);
                                }

                                var f_e_ds_values2 = e_ds_values.Select(a => new feature_info()
                                {
                                    alphabet = alphabet.name,
                                    dimension = 1,
                                    category = $@"aaindex_{aaindexes.name}",
                                    source = source.ToString(),
                                    @group = $@"aaindex_{aaindexes.name}_{sq.name}_{alphabet.name}",
                                    member = $@"{sq_index}_{alphabet_group.group_name}_{a.member_id}_{aaindex_entry.H_Accession_Number}",
                                    perspective = a.perspective_id,
                                    feature_value = a.perspective_value
                                }).ToList();

                                //if (alphabet.name == "Overall")
                                //{
                                //  Console.WriteLine();
                                //}

                                if (f_e_ds_values2.Count <= max_features) result.AddRange(f_e_ds_values2);
                            }
                        }
                    }
                }
            }

            if (calculate_aa_index_classification_data_template == null)
            {
                var template = result.Select(a => new feature_info(a) { source = "", feature_value = 0 }).ToList();
                calculate_aa_index_classification_data_template = template;
            }

            return result;
        }

        public static List<feature_info> calculate_dna_binding_prediction_data_template = null;

        public static List<feature_info> calculate_chain_dna_binding_prediction_data(List<Atom> subsequence_master_atoms, protein_data_sources source, int max_features)
        {
#if DEBUG
            if (Program.verbose_debug) Program.WriteLine($"{nameof(calculate_chain_dna_binding_prediction_data)}(List<Atom> subsequence_master_atoms, protein_data_sources source, int max_features);");
#endif

            //if (!source.Contains("protein") || !source.Contains("1d"))
            if (source != protein_data_sources.subsequence_1d)
            {
                return new List<feature_info>();
                //throw new ArgumentOutOfRangeException(nameof(source));
            }

            if (subsequence_master_atoms == null || subsequence_master_atoms.Count == 0)
            {

                if (calculate_dna_binding_prediction_data_template == null) throw new Exception();

                var template = calculate_dna_binding_prediction_data_template.Select(a => new feature_info(a)
                {
                    source = source.ToString(),
                    feature_value = 0
                }).ToList();

                return template;
            }

            var result = new List<feature_info>();

            var chain_dna_binding_prob_nr = subsequence_master_atoms.First().chain_dna_binding_prob_nr;
            var chain_dna_binding_prob_swissprot = subsequence_master_atoms.First().chain_dna_binding_prob_swissprot;
            var chain_dna_binding_prob_uniref90 = subsequence_master_atoms.First().chain_dna_binding_prob_uniref90;

            var probs = new List<(string name, double value)>()
            {
                ("nr", chain_dna_binding_prob_nr),
                ("swissprot", chain_dna_binding_prob_swissprot),
                ("uniref90", chain_dna_binding_prob_uniref90),
            };

            foreach (var prob in probs)
            {
                var f1 = new feature_info()
                {
                    alphabet = "Overall",
                    dimension = 1,
                    category = "dna_binding",
                    source = source.ToString(),
                    @group = $"dna_binding_{prob.name}",
                    member = $"{prob.name}",
                    perspective = "default",
                    feature_value = prob.value,
                };

                result.Add(f1);

                var f2 = new feature_info()
                {
                    alphabet = "Overall",
                    dimension = 1,
                    category = "dna_binding",
                    source = source.ToString(),
                    @group = $"dna_binding_all",
                    member = $"{prob.name}",
                    perspective = "default",
                    feature_value = prob.value,
                };

                result.Add(f2);
            }


            if (calculate_dna_binding_prediction_data_template == null)
            {
                var template = result.Select(a => new feature_info(a) { source = "", feature_value = 0 }).ToList();
                calculate_dna_binding_prediction_data_template = template;

            }

            return result;
        }


        public static List<feature_info> calculate_intrinsically_unordered_data_template = null;
        public static List<feature_info> calculate_intrinsically_unordered_data(List<Atom> subsequence_master_atoms, protein_data_sources source, int max_features)
        {
#if DEBUG
            if (Program.verbose_debug) Program.WriteLine($"{nameof(calculate_intrinsically_unordered_data)}(List<Atom> subsequence_master_atoms, protein_data_sources source, int max_features);");
#endif

            if (subsequence_master_atoms == null || subsequence_master_atoms.Count == 0)
            {
                if (calculate_intrinsically_unordered_data_template == null) throw new Exception();

                var template = calculate_intrinsically_unordered_data_template.Select(a => new feature_info(a)
                {
                    source = source.ToString(),
                    feature_value = 0
                }).ToList();

                return template;
            }

            var result = new List<feature_info>();

            var complete_sequence = subsequence_master_atoms;
            var split_sequence = feature_calcs.split_sequence(subsequence_master_atoms, 3, 0, false);
            var sequences = new List<(string name, List<Atom> sequence)>();
            sequences.Add(($@"unsplit", complete_sequence));
            sequences.AddRange(split_sequence.Select(a => ($@"split", a)).ToList());

            var alphabets = feature_calcs.aa_alphabets_inc_overall.ToList();
            //alphabets.Add((-1, "Overall", new List<string>() { $@"ARNDCQEGHILKMFPSTWYV" }));
            alphabets = alphabets.Where(a => !String.Equals(a.name, $@"Normal", StringComparison.InvariantCultureIgnoreCase)).ToList();

            for (var sq_index = 0; sq_index < sequences.Count; sq_index++)
            {
                var sq = sequences[sq_index];
                foreach (var alphabet in alphabets)
                {
                    var alphabet_result = new List<feature_info>();

                    foreach (var alphabet_group in alphabet.groups)
                    {
                        var iup = sq.sequence.Where(a => alphabet_group.group_amino_acids.Contains(a.amino_acid)).Select(a => a.iup_entry).ToList();

                        var short_list = iup.Select(a => a.short_type_score).ToArray();
                        var long_list = iup.Select(a => a.long_type_score).ToArray();
                        var glob_list = iup.Select(a => a.glob_type_score).ToArray();
                        var anchor2_list = iup.Select(a => a.anchor2_score).ToArray();

                        var ds_short_list = descriptive_stats.get_stat_values(short_list, $@"iup_short");
                        var ds_long_list = descriptive_stats.get_stat_values(long_list, $@"iup_long");
                        var ds_glob_list = descriptive_stats.get_stat_values(glob_list, $@"iup_glob");
                        var ds_anchor2_list = descriptive_stats.get_stat_values(anchor2_list, $@"iup_anchor2");

                        var e_ds_short_list = descriptive_stats.encode(ds_short_list);
                        var e_ds_long_list = descriptive_stats.encode(ds_long_list);
                        var e_ds_glob_list = descriptive_stats.encode(ds_glob_list);
                        var e_ds_anchor2_list = descriptive_stats.encode(ds_anchor2_list);

                        var f_e_ds_short_list = e_ds_short_list.Select(a => new feature_info()
                        {
                            alphabet = alphabet.name,
                            dimension = 1,
                            category = $@"iup",
                            source = source.ToString(),
                            @group = $@"iup_{sq.name}_short_{alphabet.name}",
                            member = $@"{sq_index}_{alphabet_group.group_name}_{a.member_id}",
                            perspective = a.perspective_id,
                            feature_value = a.perspective_value
                        }).ToList();

                        if (f_e_ds_short_list.Count <= max_features) alphabet_result.AddRange(f_e_ds_short_list);

                        var f_e_ds_long_list = e_ds_long_list.Select(a => new feature_info()
                        {
                            alphabet = alphabet.name,
                            dimension = 1,
                            category = $@"iup",
                            source = source.ToString(),
                            @group = $@"iup_{sq.name}_long_{alphabet.name}",
                            member = $@"{sq_index}_{alphabet_group.group_name}_{a.member_id}",
                            perspective = a.perspective_id,
                            feature_value = a.perspective_value
                        }).ToList();

                        if (f_e_ds_long_list.Count <= max_features) alphabet_result.AddRange(f_e_ds_long_list);

                        var f_e_ds_glob_list = e_ds_glob_list.Select(a => new feature_info()
                        {
                            alphabet = alphabet.name,
                            dimension = 1,
                            category = $@"iup",
                            source = source.ToString(),
                            @group = $@"iup_{sq.name}_glob_{alphabet.name}",
                            member = $@"{sq_index}_{alphabet_group.group_name}_{a.member_id}",
                            perspective = a.perspective_id,
                            feature_value = a.perspective_value
                        }).ToList();

                        if (f_e_ds_glob_list.Count <= max_features) alphabet_result.AddRange(f_e_ds_glob_list);

                        var f_e_ds_anchor_list = e_ds_anchor2_list.Select(a => new feature_info()
                        {
                            alphabet = alphabet.name,
                            dimension = 1,
                            category = $@"iup",
                            source = source.ToString(),
                            @group = $@"iup_{sq.name}_anchor2_{alphabet.name}",
                            member = $@"{sq_index}_{alphabet_group.group_name}_{a.member_id}",
                            perspective = a.perspective_id,
                            feature_value = a.perspective_value
                        }).ToList();

                        if (f_e_ds_anchor_list.Count <= max_features) alphabet_result.AddRange(f_e_ds_anchor_list);
                    }

                    var all = alphabet_result.Select(a => new feature_info(a) {@group = $@"iup_{sq.name}_all_{alphabet.name}",}).ToList();

                    if (all.Count <= max_features) alphabet_result.AddRange(all);


                    result.AddRange(alphabet_result);
                }
            }

            if (calculate_intrinsically_unordered_data_template == null)
            {
                var template = result.Select(a => new feature_info(a) { source = "", feature_value = 0 }).ToList();
                calculate_intrinsically_unordered_data_template = template;

            }

            return result;
        }

        public enum pssm_normalisation_methods
        {

            norm_none, norm_whole_pssm, norm_subseq, norm_encoded_parts, norm_encoded_vector
        }

        public static List<feature_info> calculate_blast_pssm_classification_data_template = null;
        public static List<feature_info> calculate_blast_pssm_classification_data(blast_pssm_options blast_pssm_options, subsequence_classification_data scd, List<Atom> subsequence_master_atoms, protein_data_sources source, int max_features)
        {
#if DEBUG
            if (Program.verbose_debug) Program.WriteLine($"{nameof(calculate_blast_pssm_classification_data)}(subsequence_classification_data scd, List<Atom> subsequence_master_atoms, protein_data_sources source, int max_features);");
#endif


            var use_databases = new List<(string, bool)>()
            {
                ("blast_pssm_nr_local_default", blast_pssm_options.db_nr_local_def),

                ("blast_pssm_nr_local_1e-4", blast_pssm_options.db_nr_local_1e_4),
                ("blast_pssm_nr_remote_1e-4", blast_pssm_options.db_nr_remote_1e_4),

                ("blast_pssm_swissprot_local_1e-4", blast_pssm_options.db_sp_local_1e_4),
                ("blast_pssm_swissprot_local_default", blast_pssm_options.db_sp_local_def),
                ("blast_pssm_swissprot_remote_1e-4", blast_pssm_options.db_sp_remote_1e_4),

                ("blast_pssm_uniref90_local_default", blast_pssm_options.db_ur90_local_def),
            };



            if (subsequence_master_atoms == null || subsequence_master_atoms.Count == 0)
            {
                if (calculate_blast_pssm_classification_data_template == null) throw new Exception();

                var template = calculate_blast_pssm_classification_data_template.Select(a => new feature_info(a)
                {
                    source = source.ToString(),
                    feature_value = 0
                }).ToList();

                return template;
            }


            var pssm_folders = Directory.GetDirectories($@"C:\betastrands_dataset\", "blast_pssm_*");
            var pssm_database_names = pssm_folders.Select(a => a.Split(new char[] { '\\', '/' }, StringSplitOptions.RemoveEmptyEntries).Last()).ToList();

            var complete_sequence = subsequence_master_atoms;


            var split_sequence = feature_calcs.split_sequence(complete_sequence, 3, 0, false);
            var sequences = new List<(string name, List<Atom> sequence)>();
            if (blast_pssm_options.make_unsplit_sequence) sequences.Add(($@"unsplit", complete_sequence));
            if (blast_pssm_options.make_split_sequence) sequences.AddRange(split_sequence.Select(a => ($@"split", a)).ToList());

            var features = new List<feature_info>();


            var ds_options_list = new List<descriptive_stats.descriptive_stats_encoding_options>();

            if (blast_pssm_options.encode_min)
            {
                ds_options_list.Add(new descriptive_stats.descriptive_stats_encoding_options(nameof(descriptive_stats.min), false) { min = true });
            }

            if (blast_pssm_options.encode_max)
            {
                ds_options_list.Add(new descriptive_stats.descriptive_stats_encoding_options(nameof(descriptive_stats.max), false) { max = true });
            }

            if (blast_pssm_options.encode_mean)
            {
                ds_options_list.Add(new descriptive_stats.descriptive_stats_encoding_options(nameof(descriptive_stats.mean_arithmetic), false) { mean_arithmetic = true });
            }

            if (blast_pssm_options.encode_mean_sd)
            {
                ds_options_list.Add(new descriptive_stats.descriptive_stats_encoding_options(nameof(descriptive_stats.mean_arithmetic), false) { mean_arithmetic = true, dev_standard = true });
            }

            if (blast_pssm_options.encode_median)
            {
                ds_options_list.Add(new descriptive_stats.descriptive_stats_encoding_options(nameof(descriptive_stats.median_q2), false) { median_q2 = true, mad_median_q2 = true });
            }

            if (blast_pssm_options.encode_mode)
            {
                ds_options_list.Add(new descriptive_stats.descriptive_stats_encoding_options(nameof(descriptive_stats.mode), false) { mode = true, mad_mode = true });
            }

            if (blast_pssm_options.encode_range)
            {
                ds_options_list.Add(new descriptive_stats.descriptive_stats_encoding_options(nameof(descriptive_stats.range), false) { range = true });
            }

            var tasks2 = new List<Task<(string database, string alphabet_name, string sequence_name, string pssm_encoding_name, List<(string member_id, string perspective_id, double perspective_value)> pssm400DT)>>();


            for (var _sq_index = 0; _sq_index < sequences.Count; _sq_index++)
            {
                var sq_index = _sq_index;

                var sq = sequences[sq_index];
                

                var all_pssm_unnormalised = sq.sequence.SelectMany(a => a.amino_acid_pssm_unnormalised).ToList();
                var all_pssm_normalised = sq.sequence.SelectMany(a => a.amino_acid_pssm_normalised).ToList();


                var alphabets = feature_calcs.aa_alphabets.ToList();

                //pssm_database_names.ForEach(a => Console.WriteLine(a));
                foreach (var _database in pssm_database_names)
                {
                    var database = _database;

                    if (!use_databases.Any(a => string.Equals(a.Item1, database, StringComparison.InvariantCultureIgnoreCase) && a.Item2)) continue;

                    var pssm_group_unnormalised = all_pssm_unnormalised.Where(a => a.database.Equals(database, StringComparison.InvariantCultureIgnoreCase)).ToList();
                    var pssm_group_normalised = all_pssm_normalised.Where(a => a.database.Equals(database, StringComparison.InvariantCultureIgnoreCase)).ToList();
                    //var pssm_entries = pssm_group.ToList();

                    var pssm_matrix_unnormalised = pssm_group_unnormalised.SelectMany(a => a.pssm_entries.Select(b => new pssm.pssm_entry(b)).ToList()).ToList();
                    var pssm_matrix_normalised = pssm_group_normalised.SelectMany(a => a.pssm_entries.Select(b => new pssm.pssm_entry(b)).ToList()).ToList();


                    // method 1: no normalisation
                    // method 2: normalise whole pssm
                    // method 3: normalise interface pssm
                    // method 4: normalise encoded pssm


                    //// !!!!!
                    //pssm_matrix_unnormalised = pssm.normalise_pssm(pssm_matrix_unnormalised);
                    //pssm_matrix_normalised = pssm.normalise_pssm(pssm_matrix_normalised);
                    /// !!!!!!

                    foreach (pssm_normalisation_methods _normalisation_method in Enum.GetValues(typeof(pssm_normalisation_methods)))
                    {
                        var normalisation_method = _normalisation_method;

                        var pssm_matrix = new List<pssm.pssm_entry>();


                        switch (normalisation_method)
                        {
                            case pssm_normalisation_methods.norm_none:
                                if (!blast_pssm_options.normalise_none) continue;
                                pssm_matrix = pssm_matrix_unnormalised;
                                break;

                            case pssm_normalisation_methods.norm_whole_pssm:
                                if (!blast_pssm_options.normalise_whole_pssm) continue;
                                pssm_matrix = pssm_matrix_normalised;
                                break;

                            case pssm_normalisation_methods.norm_subseq:
                                if (!blast_pssm_options.normalise_subsequence) continue;
                                pssm_matrix = pssm.normalise_pssm(pssm_matrix_unnormalised);
                                break;

                            case pssm_normalisation_methods.norm_encoded_parts:
                                if (!blast_pssm_options.normalise_encoded_parts) continue;
                                pssm_matrix = pssm_matrix_unnormalised;
                                break;

                            case pssm_normalisation_methods.norm_encoded_vector:
                                if (!blast_pssm_options.normalise_encoded_vector) continue;
                                pssm_matrix = pssm_matrix_unnormalised;
                                break;

                            default:
                                throw new Exception();
                        }

                        var normalisation_method_str = normalisation_method.ToString();

                        foreach (pssm.pssm_value_types _pssm_value_type in Enum.GetValues(typeof(pssm.pssm_value_types)))
                        {
                            var pssm_value_type = _pssm_value_type;

                            switch (pssm_value_type)
                            {
                                case pssm.pssm_value_types.standard:
                                    if (!blast_pssm_options.encode_standard_vector) continue;
                                    break;

                                case pssm.pssm_value_types.distances:
                                    // not supported?  need to check.
                                    continue;

                                case pssm.pssm_value_types.intervals:
                                    if (!blast_pssm_options.encoded_interval_vector) continue;
                                    break;

                                default:
                                    throw new ArgumentOutOfRangeException();
                            }

                            if (pssm_value_type == pssm.pssm_value_types.distances) continue;

                            var pssm_value_type_str = pssm_value_type.ToString();

                            foreach (var _ds_options in ds_options_list)
                            {
                                var ds_options = _ds_options;

                                var tasks = new List<Task>();

                                List<(string alphabet, List<(string col_aa, double[] values)> x)> pssm20col_values_alphabets = null;
                                List<(string alphabet, List<(string row_aa, double[] values)> x)> pssm20row_values_alphabets = null;
                                List<(string alphabet, List<(string row_aa, string col_aa, double[] values)> x)> pssm210_values_alphabets = null;
                                List<(string alphabet, List<(string row_aa, string col_aa, double[] values)> x)> pssm400_values_alphabets = null;


                                List<(string alphabet, List<(string col_aa, int lag, double[] values)> x)> pssm20colDT_values_alphabets = null;
                                List<(string alphabet, List<(string row_aa, string col_aa, int lag, double[] values)> x)> pssm210DT_values_alphabets = null;
                                List<(string alphabet, List<(string row_aa, string col_aa, int lag, double[] values)> x)> pssm400DT_values_alphabets = null;


                                var task_pssm20col_values_alphabets = !blast_pssm_options.make_standard_encoding || !blast_pssm_options.size_20 ? null : Task.Run(() => { return pssm.pssm_to_vector20col(pssm_matrix, pssm_value_type, normalisation_method == pssm_normalisation_methods.norm_encoded_parts, normalisation_method == pssm_normalisation_methods.norm_encoded_vector); });
                                var task_pssm20row_values_alphabets = !blast_pssm_options.make_standard_encoding || !blast_pssm_options.size_20 ? null : Task.Run(() => { return pssm.pssm_to_vector20row(pssm_matrix, pssm_value_type, normalisation_method == pssm_normalisation_methods.norm_encoded_parts, normalisation_method == pssm_normalisation_methods.norm_encoded_vector); });
                                var task_pssm210_values_alphabets = !blast_pssm_options.make_standard_encoding || !blast_pssm_options.size_210 ? null : Task.Run(() => { return pssm.pssm_to_vector210(pssm_matrix, pssm_value_type, normalisation_method == pssm_normalisation_methods.norm_encoded_parts, normalisation_method == pssm_normalisation_methods.norm_encoded_vector); });
                                var task_pssm400_values_alphabets = !blast_pssm_options.make_standard_encoding || !blast_pssm_options.size_400 ? null : Task.Run(() => { return pssm.pssm_to_vector400(pssm_matrix, pssm_value_type, normalisation_method == pssm_normalisation_methods.norm_encoded_parts, normalisation_method == pssm_normalisation_methods.norm_encoded_vector); });


                                if (task_pssm20col_values_alphabets != null) tasks.Add(task_pssm20col_values_alphabets);
                                if (task_pssm20row_values_alphabets != null) tasks.Add(task_pssm20row_values_alphabets);
                                if (task_pssm210_values_alphabets != null) tasks.Add(task_pssm210_values_alphabets);
                                if (task_pssm400_values_alphabets != null) tasks.Add(task_pssm400_values_alphabets);


                                var max_lag = 5;

                                var task_pssm20colDT_values_alphabets = !blast_pssm_options.make_distance_transform || !blast_pssm_options.size_20 ? null : Task.Run(() => { return pssm.pssm_to_vector20col_DT(pssm_matrix, max_lag, pssm_value_type, normalisation_method == pssm_normalisation_methods.norm_encoded_parts, normalisation_method == pssm_normalisation_methods.norm_encoded_vector); });
                                var task_pssm210DT_values_alphabets = !blast_pssm_options.make_distance_transform || !blast_pssm_options.size_210 ? null : Task.Run(() => { return pssm.pssm_to_vector210_DT(pssm_matrix, max_lag, pssm_value_type, normalisation_method == pssm_normalisation_methods.norm_encoded_parts, normalisation_method == pssm_normalisation_methods.norm_encoded_vector); });
                                var task_pssm400DT_values_alphabets = !blast_pssm_options.make_distance_transform || !blast_pssm_options.size_400 ? null : Task.Run(() => { return pssm.pssm_to_vector400_DT(pssm_matrix, max_lag, pssm_value_type, normalisation_method == pssm_normalisation_methods.norm_encoded_parts, normalisation_method == pssm_normalisation_methods.norm_encoded_vector); });

                                if (task_pssm20colDT_values_alphabets != null) tasks.Add(task_pssm20colDT_values_alphabets);
                                if (task_pssm210DT_values_alphabets != null) tasks.Add(task_pssm210DT_values_alphabets);
                                if (task_pssm400DT_values_alphabets != null) tasks.Add(task_pssm400DT_values_alphabets);

                                Task.WaitAll(tasks.ToArray<Task>());

                                pssm20col_values_alphabets = task_pssm20col_values_alphabets?.Result;
                                pssm20row_values_alphabets = task_pssm20row_values_alphabets?.Result;
                                pssm210_values_alphabets = task_pssm210_values_alphabets?.Result;
                                pssm400_values_alphabets = task_pssm400_values_alphabets?.Result;


                                pssm20colDT_values_alphabets = task_pssm20colDT_values_alphabets?.Result;
                                pssm210DT_values_alphabets = task_pssm210DT_values_alphabets?.Result;
                                pssm400DT_values_alphabets = task_pssm400DT_values_alphabets?.Result;


                                if (blast_pssm_options.make_standard_encoding)
                                {
                                    if (blast_pssm_options.size_1)
                                    {
                                        var t1 = Task.Run(() =>
                                        {
                                            var pssm1_values = pssm.pssm_to_vector1(pssm_matrix, pssm_value_type, normalisation_method == pssm_normalisation_methods.norm_encoded_vector);
                                            //if (normalise_encoding) pssm1_values = pssm.normalise_array(pssm1_values);
                                            var pssm1_ds = descriptive_stats.get_stat_values(pssm1_values, $"{sq_index}_pssm1_all_{ds_options.options_name}");

                                            var pssm1 = descriptive_stats.encode(pssm1_ds, ds_options);

                                            if (pssm1.Count > max_features) pssm1 = null;

                                            return (database, "Overall", sq.name, $"{nameof(pssm1)}_{normalisation_method_str}_{pssm_value_type_str}", pssm1);
                                        });

                                        tasks2.Add(t1);
                                    }
                                }

                                foreach (var _alphabet in alphabets)
                                {
                                    var alphabet = _alphabet;

                                    if (blast_pssm_options.make_standard_encoding)
                                    {
                                        if (blast_pssm_options.size_20)
                                        {
                                            if (pssm20col_values_alphabets != null && pssm20col_values_alphabets.Count > 0)
                                            {
                                                var pssm20col_values = pssm20col_values_alphabets.FirstOrDefault(a => a.alphabet == alphabet.name).x;

                                                if (pssm20col_values != null && pssm20col_values.Count > 0) // && pssm20col_values.Count <= max_features)
                                                {
                                                    var t = Task.Run(() =>
                                                    {
                                                        var pssm20col_ds = pssm20col_values.Select(a => descriptive_stats.get_stat_values(a.values, $"{sq_index}_pssm20_c{a.col_aa}_{ds_options.options_name}")).ToList();
                                                        var pssm20col = pssm20col_ds.SelectMany(a => a.encode(ds_options)).ToList();

                                                        if (pssm20col.Count > max_features) pssm20col = null;

                                                        return (database, alphabet.name, sq.name, $"{nameof(pssm20col)}_{normalisation_method_str}_{pssm_value_type_str}", pssm20col);
                                                    });

                                                    tasks2.Add(t);
                                                }
                                            }
                                        }

                                        if (blast_pssm_options.size_20)
                                        {
                                            if (pssm20row_values_alphabets != null && pssm20row_values_alphabets.Count > 0)
                                            {
                                                var pssm20row_values = pssm20row_values_alphabets.FirstOrDefault(a => a.alphabet == alphabet.name).x;

                                                if (pssm20row_values != null && pssm20row_values.Count > 0) // && pssm20row_values.Count <= max_features)
                                                {
                                                    var t = Task.Run(() =>
                                                    {
                                                        var pssm20row_ds = pssm20row_values.Select(a => descriptive_stats.get_stat_values(a.values, $"{sq_index}_pssm20_r{a.row_aa}_{ds_options.options_name}")).ToList();
                                                        var pssm20row = pssm20row_ds.SelectMany(a => a.encode(ds_options)).ToList();

                                                        if (pssm20row.Count > max_features) pssm20row = null;

                                                        return (database, alphabet.name, sq.name, $"{nameof(pssm20row)}_{normalisation_method_str}_{pssm_value_type_str}", pssm20row);
                                                    });

                                                    tasks2.Add(t);
                                                }
                                            }
                                        }

                                        if (blast_pssm_options.size_210)
                                        {
                                            if (pssm210_values_alphabets != null && pssm210_values_alphabets.Count > 0)
                                            {
                                                var pssm210_values = pssm210_values_alphabets.FirstOrDefault(a => a.alphabet == alphabet.name).x;

                                                if (pssm210_values != null && pssm210_values.Count > 0) // && pssm210_values.Count <= max_features)
                                                {
                                                    var t = Task.Run(() =>
                                                    {
                                                        var pssm210_ds = pssm210_values.Select(a => descriptive_stats.get_stat_values(a.values, $"{sq_index}_pssm210_c{a.col_aa}_r{a.row_aa}_{ds_options.options_name}")).ToList();
                                                        var pssm210 = pssm210_ds.SelectMany(a => a.encode(ds_options)).ToList();

                                                        if (pssm210.Count > max_features) pssm210 = null;

                                                        return (database, alphabet.name, sq.name, $"{nameof(pssm210)}_{normalisation_method_str}_{pssm_value_type_str}", pssm210);
                                                    });

                                                    tasks2.Add(t);
                                                }
                                            }
                                        }

                                        if (blast_pssm_options.size_400)
                                        {
                                            if (pssm400_values_alphabets != null && pssm400_values_alphabets.Count > 0)
                                            {
                                                var pssm400_values = pssm400_values_alphabets.FirstOrDefault(a => a.alphabet == alphabet.name).x;

                                                if (pssm400_values != null && pssm400_values.Count > 0) // && pssm400_values.Count <= max_features)
                                                {
                                                    var t = Task.Run(() =>
                                                    {
                                                        var pssm400_ds = pssm400_values.Select(a => descriptive_stats.get_stat_values(a.values, $"{sq_index}_pssm400_c{a.col_aa}_r{a.row_aa}_{ds_options.options_name}")).ToList();
                                                        var pssm400 = pssm400_ds.SelectMany(a => a.encode(ds_options)).ToList();

                                                        if (pssm400.Count > max_features) pssm400 = null;

                                                        return (database, alphabet.name, sq.name, $"{nameof(pssm400)}_{normalisation_method_str}_{pssm_value_type_str}", pssm400);
                                                    });

                                                    tasks2.Add(t);
                                                }
                                            }
                                        }
                                    }


                                    if (blast_pssm_options.make_distance_transform)
                                    {
                                        if (blast_pssm_options.size_20)
                                        {
                                            if (pssm20colDT_values_alphabets != null && pssm20colDT_values_alphabets.Count > 0)
                                            {
                                                var pssm20colDT_values_groups = pssm20colDT_values_alphabets.FirstOrDefault(a => a.alphabet == alphabet.name).x.GroupBy(a => a.lag).ToList();

                                                if (pssm20colDT_values_groups != null && pssm20colDT_values_groups.Count > 0)
                                                {
                                                    foreach (var _pssm20colDT_values_group in pssm20colDT_values_groups)
                                                    {
                                                        var pssm20colDT_values_group = _pssm20colDT_values_group;

                                                        var pssm20colDT_lag = pssm20colDT_values_group.Key;
                                                        var pssm20colDT_values = pssm20colDT_values_group.ToList();

                                                        if (pssm20colDT_values.Count > 0) // && pssm20colDT_values.Count <= max_features)
                                                        {
                                                            var t = Task.Run(() =>
                                                            {
                                                                var pssm20colDT_ds = pssm20colDT_values.Select(a => descriptive_stats.get_stat_values(a.values, $"{sq_index}_pssm20colDT_lag{pssm20colDT_lag}_c{a.col_aa}_rx_{ds_options.options_name}")).ToList();
                                                                var pssm20colDT = pssm20colDT_ds.SelectMany(a => a.encode(ds_options)).ToList();

                                                                if (pssm20colDT.Count > max_features) pssm20colDT = null;

                                                                return (database, alphabet.name, sq.name, $"{nameof(pssm20colDT)}_lag{pssm20colDT_lag}_{normalisation_method_str}_{pssm_value_type_str}", pssm20colDT);
                                                            });
                                                            tasks2.Add(t);
                                                        }
                                                    }
                                                }
                                            }
                                        }

                                        if (blast_pssm_options.size_210)
                                        {
                                            if (pssm210DT_values_alphabets != null && pssm210DT_values_alphabets.Count > 0)
                                            {
                                                var pssm210DT_values_groups = pssm210DT_values_alphabets.FirstOrDefault(a => a.alphabet == alphabet.name).x.GroupBy(a => a.lag).ToList();

                                                if (pssm210DT_values_groups != null && pssm210DT_values_groups.Count > 0)
                                                {
                                                    foreach (var _pssm210DT_values_group in pssm210DT_values_groups)
                                                    {
                                                        var pssm210DT_values_group = _pssm210DT_values_group;

                                                        var pssm210DT_lag = pssm210DT_values_group.Key;
                                                        var pssm210DT_values = pssm210DT_values_group.ToList();
                                                        if (pssm210DT_values.Count > 0) // && pssm210DT_values.Count <= max_features)
                                                        {
                                                            var t = Task.Run(() =>
                                                            {
                                                                var pssm210DT_ds = pssm210DT_values.Select(a => descriptive_stats.get_stat_values(a.values, $"{sq_index}_pssm210DT_lag{pssm210DT_lag}_c{a.col_aa}_r{a.row_aa}_{ds_options.options_name}")).ToList();
                                                                var pssm210DT = pssm210DT_ds.SelectMany(a => a.encode(ds_options)).ToList();

                                                                if (pssm210DT.Count > max_features) pssm210DT = null;

                                                                return (database, alphabet.name, sq.name, $"{nameof(pssm210DT)}_lag{pssm210DT_lag}_{normalisation_method_str}_{pssm_value_type_str}", pssm210DT);
                                                            });

                                                            tasks2.Add(t);
                                                        }
                                                    }
                                                }
                                            }
                                        }

                                        if (blast_pssm_options.size_400)
                                        {
                                            if (pssm400DT_values_alphabets != null && pssm400DT_values_alphabets.Count > 0)
                                            {
                                                var pssm400DT_values_groups = pssm400DT_values_alphabets.FirstOrDefault(a => a.alphabet == alphabet.name).x.GroupBy(a => a.lag).ToList();

                                                if (pssm400DT_values_groups != null && pssm400DT_values_groups.Count > 0)
                                                {
                                                    foreach (var _pssm400DT_values_group in pssm400DT_values_groups)
                                                    {
                                                        var pssm400DT_values_group = _pssm400DT_values_group;
                                                        var pssm400DT_lag = pssm400DT_values_group.Key;
                                                        var pssm400DT_values = pssm400DT_values_group.ToList();

                                                        if (pssm400DT_values.Count > 0) // && pssm400DT_values.Count <= max_features)
                                                        {
                                                            var t = Task.Run(() =>
                                                            {
                                                                var pssm400DT_ds = pssm400DT_values.Select(a => descriptive_stats.get_stat_values(a.values, $"{sq_index}_pssm400DT_lag{pssm400DT_lag}_c{a.col_aa}_r{a.row_aa}_{ds_options.options_name}")).ToList();
                                                                var pssm400DT = pssm400DT_ds.SelectMany(a => a.encode(ds_options)).ToList();

                                                                if (pssm400DT.Count > max_features) pssm400DT = null;

                                                                return (database, alphabet.name, sq.name, $"{nameof(pssm400DT)}_lag{pssm400DT_lag}_{normalisation_method_str}_{pssm_value_type_str}", pssm400DT);
                                                            });

                                                            tasks2.Add(t);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            var pssm_list = new List<(string database, string alphabet, string sequence_name, string pssm_encoding_name, List<(string member_id, string perspective_id, double perspective_value)> values)>();

            if (tasks2 != null && tasks2.Count > 0)
            {
                Task.WaitAll(tasks2.ToArray<Task>());

                var result = tasks2.Select(a => a.Result).ToList();

                pssm_list.AddRange(result);
            }

            foreach ((string database, string alphabet, string sequence_name, string pssm_encoding_name, List<(string member_id, string perspective_id, double perspective_value)> values) x in pssm_list)
            {
                if (x.values == null || x.values.Count == 0) continue;

                if (x.values.Count > max_features) continue;
                
                foreach ((string member_id, string perspective_id, double perspective_value) y in x.values)
                {
                    features.Add(new feature_info()
                    {
                        alphabet = x.alphabet,
                        dimension = 1,
                        category = $"blast_pssm",
                        source = source.ToString(),
                        group = $"{x.database}_{x.sequence_name}_{x.alphabet}_{x.pssm_encoding_name}",
                        member = y.member_id,
                        perspective = y.perspective_id,
                        feature_value = y.perspective_value
                    });
                }
            }

            if (calculate_blast_pssm_classification_data_template == null)
            {
                var template = features.Select(a => new feature_info(a) { source = "", feature_value = 0 }).ToList();
                calculate_blast_pssm_classification_data_template = template;
            }

            
            //var g = features.Select(a => a.@group).Distinct().Count();

            //Console.WriteLine("blast pssm features: " + features.Count + " groups: " + g);
            return features;
        }

        public static List<feature_info> calculate_sasa_classification_data_template = null;

        public static List<feature_info> calculate_sasa_classification_data(List<Atom> subsequence_master_atoms, protein_data_sources source, int max_features)
        {
#if DEBUG
            if (Program.verbose_debug) Program.WriteLine($"{nameof(calculate_sasa_classification_data)}(List<Atom> subsequence_master_atoms, protein_data_sources source, int max_features);");
#endif

            if (subsequence_master_atoms == null || subsequence_master_atoms.Count == 0)
            {
                if (calculate_sasa_classification_data_template == null) throw new Exception();

                var template = calculate_sasa_classification_data_template.Select(a => new feature_info(a)
                {
                    source = source.ToString(),
                    feature_value = 0
                }).ToList();

                return template;
            }

            var features = new List<feature_info>();

            var complete_sequence = subsequence_master_atoms;
            var split_sequence = feature_calcs.split_sequence(complete_sequence, 3, 0, false);
            var sequences = new List<(string name, List<Atom> sequence)>();
            sequences.Add(($@"unsplit", complete_sequence));
            sequences.AddRange(split_sequence.Select(a => ($@"split", a)).ToList());


            var alphabets = feature_calcs.aa_alphabets_inc_overall.ToList();
            //alphabets.Add((-1, "Overall", new List<string>() { "ARNDCQEGHILKMFPSTWYV" }));
            //alphabets = alphabets.Where(a => !String.Equals(a.name, "Normal", StringComparison.InvariantCultureIgnoreCase)).ToList();

            for (var _sq_index = 0; _sq_index < sequences.Count; _sq_index++)
            {
                var sq_index = _sq_index;
                var sq = sequences[sq_index];

                foreach (var alphabet in alphabets)
                {
                    var all_rel = new List<feature_info>(); // rel & (s | l)
                    var all_abs = new List<feature_info>(); // abs & (s | l)
                    var all_algorithm_s = new List<feature_info>(); // S & (abs | rel)
                    var all_algorithm_l = new List<feature_info>(); // L & (abs | rel)
                    var all_rel_abs_s_l = new List<feature_info>();

                    foreach (var alphabet_group in alphabet.groups)
                    {
                        var sequence_sasa_values = sq.sequence.Where(a => alphabet_group.group_amino_acids.Contains(a.amino_acid)).Select(a => (L_all_atoms_abs: a.RSA_L.all_atoms_abs, L_all_polar_abs: a.RSA_L.all_polar_abs, L_main_chain_abs: a.RSA_L.main_chain_abs, L_total_side_abs: a.RSA_L.total_side_abs, L_non_polar_abs: a.RSA_L.non_polar_abs, L_all_atoms_rel: a.RSA_L.all_atoms_rel, L_all_polar_rel: a.RSA_L.all_polar_rel, L_main_chain_rel: a.RSA_L.main_chain_rel, L_total_side_rel: a.RSA_L.total_side_rel, L_non_polar_rel: a.RSA_L.non_polar_rel, S_all_atoms_abs: a.RSA_S.all_atoms_abs, S_all_polar_abs: a.RSA_S.all_polar_abs, S_main_chain_abs: a.RSA_S.main_chain_abs, S_total_side_abs: a.RSA_S.total_side_abs, S_non_polar_abs: a.RSA_S.non_polar_abs, S_all_atoms_rel: a.RSA_S.all_atoms_rel, S_all_polar_rel: a.RSA_S.all_polar_rel, S_main_chain_rel: a.RSA_S.main_chain_rel, S_total_side_rel: a.RSA_S.total_side_rel, S_non_polar_rel: a.RSA_S.non_polar_rel)).ToList();

                        var all_atoms_abs_L = sequence_sasa_values.Select(a => a.L_all_atoms_abs).ToArray();
                        var all_polar_abs_L = sequence_sasa_values.Select(a => a.L_all_polar_abs).ToArray();
                        var main_chain_abs_L = sequence_sasa_values.Select(a => a.L_main_chain_abs).ToArray();
                        var total_side_abs_L = sequence_sasa_values.Select(a => a.L_total_side_abs).ToArray();
                        var non_polar_abs_L = sequence_sasa_values.Select(a => a.L_non_polar_abs).ToArray();

                        var all_atoms_rel_L = sequence_sasa_values.Select(a => a.L_all_atoms_rel).ToArray();
                        var all_polar_rel_L = sequence_sasa_values.Select(a => a.L_all_polar_rel).ToArray();
                        var main_chain_rel_L = sequence_sasa_values.Select(a => a.L_main_chain_rel).ToArray();
                        var total_side_rel_L = sequence_sasa_values.Select(a => a.L_total_side_rel).ToArray();
                        var non_polar_rel_L = sequence_sasa_values.Select(a => a.L_non_polar_rel).ToArray();

                        var all_atoms_abs_S = sequence_sasa_values.Select(a => a.S_all_atoms_abs).ToArray();
                        var all_polar_abs_S = sequence_sasa_values.Select(a => a.S_all_polar_abs).ToArray();
                        var main_chain_abs_S = sequence_sasa_values.Select(a => a.S_main_chain_abs).ToArray();
                        var total_side_abs_S = sequence_sasa_values.Select(a => a.S_total_side_abs).ToArray();
                        var non_polar_abs_S = sequence_sasa_values.Select(a => a.S_non_polar_abs).ToArray();

                        var all_atoms_rel_S = sequence_sasa_values.Select(a => a.S_all_atoms_rel).ToArray();
                        var all_polar_rel_S = sequence_sasa_values.Select(a => a.S_all_polar_rel).ToArray();
                        var main_chain_rel_S = sequence_sasa_values.Select(a => a.S_main_chain_rel).ToArray();
                        var total_side_rel_S = sequence_sasa_values.Select(a => a.S_total_side_rel).ToArray();
                        var non_polar_rel_S = sequence_sasa_values.Select(a => a.S_non_polar_rel).ToArray();

                        var all = new List<(string algo, string abs_or_rel, string sasa_type, double[] values, descriptive_stats ds_values)>();


                        all.Add(("L", "abs", "all_atoms", all_atoms_abs_L, descriptive_stats.get_stat_values(all_atoms_abs_L, nameof(all_atoms_abs_L))));
                        all.Add(("L", "abs", "all_polar", all_polar_abs_L, descriptive_stats.get_stat_values(all_polar_abs_L, nameof(all_polar_abs_L))));
                        all.Add(("L", "abs", "main_chain", main_chain_abs_L, descriptive_stats.get_stat_values(main_chain_abs_L, nameof(main_chain_abs_L))));
                        all.Add(("L", "abs", "total_side", total_side_abs_L, descriptive_stats.get_stat_values(total_side_abs_L, nameof(total_side_abs_L))));
                        all.Add(("L", "abs", "non_polar", non_polar_abs_L, descriptive_stats.get_stat_values(non_polar_abs_L, nameof(non_polar_abs_L))));

                        all.Add(("L", "rel", "all_atoms", all_atoms_rel_L, descriptive_stats.get_stat_values(all_atoms_rel_L, nameof(all_atoms_rel_L))));
                        all.Add(("L", "rel", "all_polar", all_polar_rel_L, descriptive_stats.get_stat_values(all_polar_rel_L, nameof(all_polar_rel_L))));
                        all.Add(("L", "rel", "main_chain", main_chain_rel_L, descriptive_stats.get_stat_values(main_chain_rel_L, nameof(main_chain_rel_L))));
                        all.Add(("L", "rel", "total_side", total_side_rel_L, descriptive_stats.get_stat_values(total_side_rel_L, nameof(total_side_rel_L))));
                        all.Add(("L", "rel", "non_polar", non_polar_rel_L, descriptive_stats.get_stat_values(non_polar_rel_L, nameof(non_polar_rel_L))));

                        all.Add(("S", "abs", "all_atoms", all_atoms_abs_S, descriptive_stats.get_stat_values(all_atoms_abs_S, nameof(all_atoms_abs_S))));
                        all.Add(("S", "abs", "all_polar", all_polar_abs_S, descriptive_stats.get_stat_values(all_polar_abs_S, nameof(all_polar_abs_S))));
                        all.Add(("S", "abs", "main_chain", main_chain_abs_S, descriptive_stats.get_stat_values(main_chain_abs_S, nameof(main_chain_abs_S))));
                        all.Add(("S", "abs", "total_side", total_side_abs_S, descriptive_stats.get_stat_values(total_side_abs_S, nameof(total_side_abs_S))));
                        all.Add(("S", "abs", "non_polar", non_polar_abs_S, descriptive_stats.get_stat_values(non_polar_abs_S, nameof(non_polar_abs_S))));

                        all.Add(("S", "rel", "all_atoms", all_atoms_rel_S, descriptive_stats.get_stat_values(all_atoms_rel_S, nameof(all_atoms_rel_S))));
                        all.Add(("S", "rel", "all_polar", all_polar_rel_S, descriptive_stats.get_stat_values(all_polar_rel_S, nameof(all_polar_rel_S))));
                        all.Add(("S", "rel", "main_chain", main_chain_rel_S, descriptive_stats.get_stat_values(main_chain_rel_S, nameof(main_chain_rel_S))));
                        all.Add(("S", "rel", "total_side", total_side_rel_S, descriptive_stats.get_stat_values(total_side_rel_S, nameof(total_side_rel_S))));
                        all.Add(("S", "rel", "non_polar", non_polar_rel_S, descriptive_stats.get_stat_values(non_polar_rel_S, nameof(non_polar_rel_S))));

                        foreach (var x in all.GroupBy(a => (a.algo, a.abs_or_rel)).ToList())
                        {
                            var algo = x.Key.algo;
                            var abs_or_rel = x.Key.abs_or_rel;

                            var group_list = x.ToList(); // e.g. list of all 'S' & 'rel

                            var e = group_list.Select(a => descriptive_stats.encode(a.ds_values)).ToList();

                            var f = e.SelectMany(a => a.Select(b => new feature_info()
                            {
                                alphabet = alphabet.name,
                                dimension = 3,
                                category = $@"sasa",
                                source = source.ToString(),
                                @group = $@"sasa_{sq.name}_{algo}_{abs_or_rel}_{alphabet.name}",
                                member = $@"{sq_index}_{b.member_id}_{alphabet_group}",
                                perspective = b.perspective_id,
                                feature_value = b.perspective_value
                            })).ToList();

                            if (f.Count <= max_features) features.AddRange(f);

                            if (abs_or_rel == $@"abs") all_abs.AddRange(f.Select(a => new feature_info(a) {@group = $@"sasa_{sq.name}_all_{abs_or_rel}_{alphabet.name}"}).ToList());
                            if (abs_or_rel == $@"rel") all_rel.AddRange(f.Select(a => new feature_info(a) {@group = $@"sasa_{sq.name}_all_{abs_or_rel}_{alphabet.name}"}).ToList());
                            if (algo == $@"L") all_algorithm_l.AddRange(f.Select(a => new feature_info(a) {@group = $@"sasa_{sq.name}_all_{algo}_{alphabet.name}"}).ToList());
                            if (algo == $@"S") all_algorithm_s.AddRange(f.Select(a => new feature_info(a) {@group = $@"sasa_{sq.name}_all_{algo}_{alphabet.name}"}).ToList());
                            all_rel_abs_s_l.AddRange(f.Select(a => new feature_info(a) {@group = $@"sasa_{sq.name}_all_{alphabet.name}"}).ToList());
                        }
                    }

                    if (all_abs.Count <= max_features) features.AddRange(all_abs);
                    if (all_rel.Count <= max_features) features.AddRange(all_rel);
                    if (all_algorithm_l.Count <= max_features) features.AddRange(all_algorithm_l);
                    if (all_algorithm_s.Count <= max_features) features.AddRange(all_algorithm_s);
                    if (all_rel_abs_s_l.Count <= max_features) features.AddRange(all_rel_abs_s_l);
                }
            }

            if (calculate_sasa_classification_data_template == null)
            {
                var template = features.Select(a => new feature_info(a) { source = "", feature_value = 0 }).ToList();
                calculate_sasa_classification_data_template = template;
            }

            return features;
        }

        public static List<feature_info> calculate_tortuosity_classification_data_template = null;

        public static List<feature_info> calculate_tortuosity_classification_data(subsequence_classification_data scd, List<Atom> subsequence_master_atoms, protein_data_sources source, int max_features)
        {
#if DEBUG
            if (Program.verbose_debug) Program.WriteLine($"{nameof(calculate_tortuosity_classification_data)}(subsequence_classification_data scd, List<Atom> subsequence_master_atoms, protein_data_sources source, int max_features);");
#endif

            if (subsequence_master_atoms == null || subsequence_master_atoms.Count == 0)
            {
                if (calculate_tortuosity_classification_data_template == null) throw new Exception();

                var template = calculate_tortuosity_classification_data_template.Select(a => new feature_info(a)
                {
                    source = source.ToString(),
                    feature_value = 0
                }).ToList();

                return template;
            }

            var features = new List<feature_info>();

            var complete_sequence = subsequence_master_atoms;
            var split_sequence = feature_calcs.split_sequence(complete_sequence, 3, 0, false);
            var sequences = new List<(string name, List<Atom> sequence)>();
            sequences.Add(($@"unsplit", complete_sequence));
            sequences.AddRange(split_sequence.Select(a => ($@"split", a)).ToList());


            var pdb_unsplit_features = new List<feature_info>();


            for (var sq_index = 0; sq_index < sequences.Count; sq_index++)
            {
                var sq = sequences[sq_index];
                var feats = new List<feature_info>();

                var tortuosity1 = Atom.measure_tortuosity1(sq.sequence);

                var tortuosity2 = Atom.measure_tortuosity2(sq.sequence);


                // tortuosity 1

                var x0 = new feature_info()
                {
                    alphabet = $@"Overall",
                    dimension = 3,
                    category = $@"geometry",
                    source = source.ToString(),
                    @group = $@"geometry_{sq.name}_tortuosity1",
                    member = $@"{sq_index}_default",
                    perspective = $@"default",
                    feature_value = tortuosity1.tortuosity1
                };
                feats.Add(x0);

                // tortuosity 2
                var tortuosity2_tortuosity_stat_values_encoded = descriptive_stats.encode(tortuosity2.tortuosity_stat_values);
                var x1 = tortuosity2_tortuosity_stat_values_encoded.Select(a => new feature_info()
                {
                    alphabet = $@"Overall",
                    dimension = 3,
                    category = $@"geometry",
                    source = source.ToString(),
                    @group = $@"geometry_{sq.name}_tortuosity2",
                    member = $@"{sq_index}_{a.member_id}",
                    perspective = a.perspective_id,
                    feature_value = a.perspective_value
                }).ToList();
                feats.AddRange(x1);

                // tort 1 and tort 2

                var x2 = new feature_info()
                {
                    alphabet = $@"Overall",
                    dimension = 3,
                    category = $@"geometry",
                    source = source.ToString(),
                    @group = $@"geometry_{sq.name}_tortuosity1_and_tortuosity2",
                    member = $@"{sq_index}_{nameof(tortuosity1)}",
                    perspective = $@"default",
                    feature_value = tortuosity1.tortuosity1
                };
                feats.Add(x2);

                var x3 = tortuosity2_tortuosity_stat_values_encoded.Select(a => new feature_info()
                {
                    alphabet = $@"Overall",
                    dimension = 3,
                    category = $@"geometry",
                    source = source.ToString(),
                    @group = $@"geometry_{sq.name}_tortuosity1_and_tortuosity2",
                    member = $@"{sq_index}_{a.member_id}",
                    perspective = a.perspective_id,
                    feature_value = a.perspective_value
                }).ToList();
                feats.AddRange(x3);

                // displacement length (global)
                var x4 = new feature_info()
                {
                    alphabet = $@"Overall",
                    dimension = 3,
                    category = $@"geometry",
                    source = source.ToString(),
                    @group = $@"geometry_{sq.name}_displacement_3d_global",
                    member = $@"{sq_index}_default",
                    perspective = $@"default",
                    feature_value = tortuosity1.displacement
                };
                feats.Add(x4);

                // curve length (global)
                var x5 = new feature_info()
                {
                    alphabet = $@"Overall",
                    dimension = 3,
                    category = $@"geometry",
                    source = source.ToString(),
                    @group = $@"geometry_{sq.name}_peptide_length_3d_global",
                    member = $@"{sq_index}_default",
                    perspective = $@"default",
                    feature_value = tortuosity1.distance_of_curve
                };
                feats.Add(x5);

                // average displacement length (local)
                var tortuosity2_displacements_encoded = descriptive_stats.encode(tortuosity2.displacement_stat_values);
                var x6 = tortuosity2_displacements_encoded.Select(a => new feature_info()
                {
                    alphabet = $@"Overall",
                    dimension = 3,
                    category = $@"geometry",
                    source = source.ToString(),
                    @group = $@"geometry_{sq.name}_displacement_3d_local",
                    member = $@"{sq_index}_{a.member_id}",
                    perspective = a.perspective_id,
                    feature_value = a.perspective_value
                }).ToList();
                feats.AddRange(x6);

                // average curve length (local)
                var tortuosity2_curves_encoded = descriptive_stats.encode(tortuosity2.curve_stat_values);
                var x7 = tortuosity2_curves_encoded.Select(a => new feature_info()
                {
                    alphabet = $@"Overall",
                    dimension = 3,
                    category = $@"geometry",
                    source = source.ToString(),
                    @group = $@"geometry_{sq.name}_peptide_length_3d_local",
                    member = $@"{sq_index}_{a.member_id}",
                    perspective = a.perspective_id,
                    feature_value = a.perspective_value
                }).ToList();
                feats.AddRange(x7);

                if (sq.name.StartsWith("pdb_unsplit"))
                {
                    pdb_unsplit_features = feats;
                }

                else if (sq.name.StartsWith("unsplit"))
                {
                    var rel_feats = feats.Select((a, i) => new feature_info(a)
                    {
                        @group = a.@group.Replace($"geometry_{sq.name}", $"geometry_rel_{sq.name}"),
                        feature_value = pdb_unsplit_features[i].feature_value != 0 ? a.feature_value / pdb_unsplit_features[i].feature_value : 0
                    }).ToList();

                    feats.AddRange(rel_feats);
                }

                features.AddRange(feats);
            }

            if (calculate_tortuosity_classification_data_template == null)
            {
                var template = features.Select(a => new feature_info(a) { source = "", feature_value = 0 }).ToList();
                calculate_tortuosity_classification_data_template = template;
            }


            return features;
        }

        //calculate_aa_aa_distances_classification_data

        public static List<feature_info> calculate_aa_aa_distances_classification_data_template = null;

        public static List<feature_info> calculate_aa_aa_distances_classification_data(List<Atom> subsequence_master_atoms, protein_data_sources source, int max_features)
        {
#if DEBUG
            if (Program.verbose_debug) Program.WriteLine($"{nameof(calculate_aa_aa_distances_classification_data)}List<Atom> subsequence_master_atoms, protein_data_sources source, int max_features);");
#endif

            if (subsequence_master_atoms == null || subsequence_master_atoms.Count == 0)
            {
                if (calculate_aa_aa_distances_classification_data_template == null) throw new Exception();

                var template = calculate_aa_aa_distances_classification_data_template.Select(a => new feature_info(a) {source = source.ToString(), feature_value = 0}).ToList();

                return template;
            }

            var features = new List<feature_info>();

            var distances = Atom.get_master_contacts(null, subsequence_master_atoms);
            var aa_distances = distances.Select(a => (aa1: a.atom1.amino_acid, aa2: a.atom2.amino_acid, distance: a.distance)).ToList();

            foreach (var alphabet in feature_calcs.aa_alphabets)
            {

                var fa = new List<feature_info>();

                for (var group_index1 = 0; group_index1 < alphabet.groups.Count; group_index1++)
                {
                    var group1 = alphabet.groups[group_index1];

                    for (var group_index2 = 0; group_index2 < alphabet.groups.Count; group_index2++)
                    {
                        if (group_index2 > group_index1) continue;
                        

                        var group2 = alphabet.groups[group_index2];

                        var dist_list = aa_distances.Where(a => (group1.group_amino_acids.Contains(a.aa1) && group2.group_amino_acids.Contains(a.aa2)) || (group1.group_amino_acids.Contains(a.aa2) && group2.group_amino_acids.Contains(a.aa1))).Select(a => a.distance).ToArray();

                        var ds = descriptive_stats.get_stat_values(dist_list, $"{group1.group_name}_{group2.group_name}");
                        var dse = descriptive_stats.encode(ds);

                        var f = dse.Select(a => new feature_info()
                        {
                            alphabet = alphabet.name,
                            dimension = 3,
                            category = "aa_to_aa_distances",
                            source = source.ToString(),
                            @group = $@"aa_to_aa_distances_{alphabet.name}",
                            member = a.member_id,
                            perspective = a.perspective_id,
                            feature_value = a.perspective_value,
                        }).ToList();

                        fa.AddRange(f);
                    }
                }

                if (fa.Count <= max_features) features.AddRange(fa);
            }

            if (calculate_aa_aa_distances_classification_data_template == null)
            {
                var template = features.Select(a => new feature_info(a) { source = "", feature_value = 0 }).ToList();
                calculate_aa_aa_distances_classification_data_template = template;
            }

            return features;

        }


        public static List<feature_info> calculate_intramolecular_classification_data_template = null;
        public static List<feature_info> calculate_intramolecular_classification_data(List<Atom> subsequence_master_atoms, protein_data_sources source, int max_features)
        {
#if DEBUG
            if (Program.verbose_debug) Program.WriteLine($"{nameof(calculate_intramolecular_classification_data)}List<Atom> subsequence_master_atoms, protein_data_sources source, int max_features);");
#endif

            if (subsequence_master_atoms == null || subsequence_master_atoms.Count == 0)
            {
                if (calculate_intramolecular_classification_data_template == null) throw new Exception();

                var template = calculate_intramolecular_classification_data_template.Select(a => new feature_info(a)
                {
                    source = source.ToString(),
                    feature_value = 0
                }).ToList();

                return template;
            }

            var make_contact_distance_features = true;
            var make_contact_count_features = true;

            var features = new List<feature_info>();

            // average number of intra-molecular contacts per atom
            var x = subsequence_master_atoms.Select(a => (double)(a.contact_map_intramolecular?.Count ?? 0)).ToArray();
            var intramolecular_contact_count = descriptive_stats.get_stat_values(x, "intramolecular_contact_count");

            // average distance between intra-molecular contacts
            var y = subsequence_master_atoms.Where(a => a.contact_map_intramolecular != null && a.contact_map_intramolecular.Count > 0).SelectMany(a => a.contact_map_intramolecular.Select(b => b.distance).ToList()).ToArray();
            var intramolecular_contact_distance = descriptive_stats.get_stat_values(y, "intramolecular_contact_distance");

            if (make_contact_distance_features)
            {
                var intramolecular_contact_distance_encoded = descriptive_stats.encode(intramolecular_contact_distance);

                var x0 = intramolecular_contact_distance_encoded.Select(a => new feature_info()
                {
                    alphabet = "Overall",
                    dimension = 3,
                    category = "geometry",

                    source = source.ToString(),
                    group = "geometry_" + nameof(intramolecular_contact_distance),
                    member = a.member_id,
                    perspective = a.perspective_id,
                    feature_value = a.perspective_value
                }).ToList();
                features.AddRange(x0);

                if (make_contact_distance_features && make_contact_count_features)
                {
                    var x1 = intramolecular_contact_distance_encoded.Select(a => new feature_info()
                    {
                        alphabet = "Overall",
                        dimension = 3,
                        category = "geometry",

                        source = source.ToString(),
                        group = "geometry_intramolecular_contact_count_and_distance",
                        member = a.member_id,
                        perspective = a.perspective_id,
                        feature_value = a.perspective_value
                    }).ToList();
                    features.AddRange(x1);
                }
            }

            if (make_contact_count_features)
            {
                var intramolecular_contact_count_encoded = descriptive_stats.encode(intramolecular_contact_count);
                var x0 = intramolecular_contact_count_encoded.Select(a => new feature_info()
                {
                    alphabet = "Overall",
                    dimension = 3,
                    category = "geometry",

                    source = source.ToString(),
                    group = "geometry_" + nameof(intramolecular_contact_count),
                    member = a.member_id,
                    perspective = a.perspective_id,
                    feature_value = a.perspective_value
                }).ToList();
                features.AddRange(x0);


                if (make_contact_distance_features && make_contact_count_features)
                {
                    var x1 = intramolecular_contact_count_encoded.Select(a => new feature_info()
                    {
                        alphabet = "Overall",
                        dimension = 3,
                        category = "geometry",

                        source = source.ToString(),
                        group = "geometry_intramolecular_contact_count_and_distance",
                        member = a.member_id,
                        perspective = a.perspective_id,
                        feature_value = a.perspective_value
                    }).ToList();

                    features.AddRange(x1);
                }
            }

            if (calculate_intramolecular_classification_data_template == null)
            {
                var template = features.Select(a => new feature_info(a) { source = "", feature_value = 0 }).ToList();
                calculate_intramolecular_classification_data_template = template;
            }

            return features;
        }


        public static int total_pse_aac_sequence_classification_data = -1;
        public static int total_sable_sequence_classification_data = -1;
        public static int total_mpsa_classification_data = -1;
        public static int total_blast_pssm_subsequence_classification_data = -1;
        //public static int total_blast_pssm_protein_classification_data = -1;
        public static int total_aa_index_classification_data = -1;
        public static int total_sequence_geometry_classification_data = -1;
        public static int total_intrinsically_unordered_data = -1;
        public static int total_dna_binding_prediction_data = -1;
        public static int total_r_peptides_prediction_data = -1;

        public static int total_subsequence_1d_classification_data = -1;
        public static int total_neighbourhood_1d_classification_data = -1;
        public static int total_protein_1d_classification_data = -1;

        //public static int total_protein_pse_ssc_dssp_classification_data = -1;
        public static int total_pse_ssc_dssp_classification_data = -1;
        public static int total_protein_dssp_dist_classification_data = -1;

        public static int total_foldx_classification_subsequence_3d_data = -1;
        public static int total_foldx_classification_neighbourhood_3d_data = -1;
        public static int total_foldx_classification_protein_3d_data = -1;

        public static int total_ring_classification_data = -1;
        public static int total_sasa_classification_data = -1;
        public static int total_tortuosity_classification_data = -1;
        public static int total_intramolecular_classification_data = -1;
        public static int total_aa_aa_distances_classification_data = -1;

        public static int total_subsequence_3d_classification_data = -1;
        public static int total_neighbourhood_3d_classification_data = -1;
        public static int total_protein_3d_classification_data = -1;

        public class feature_types_1d
        {
            public bool pse_aac_sequence_classification_data = false;
            public bool sable_classification_data = false;
            public bool mpsa_classification_data_subsequence = false;
            //public bool mpsa_classification_data_protein = false;
            public bool blast_pssm_subsequence_classification_data = false;
            public bool aa_index_classification_data = false;
            public bool sequence_geometry_classification_data = false;
            public bool intrinsically_unordered_data = false;
            public bool dna_binding_prediction_data = false;
            public bool r_peptides = false;

            public feature_types_1d()
            {

            }
            public feature_types_1d(bool enable)
            {
                pse_aac_sequence_classification_data = enable;
                sable_classification_data = enable;
                mpsa_classification_data_subsequence = enable;
                //mpsa_classification_data_protein = enable;
                blast_pssm_subsequence_classification_data = enable;
                aa_index_classification_data = enable;
                sequence_geometry_classification_data = enable;
                intrinsically_unordered_data = enable;
                dna_binding_prediction_data = enable;
                r_peptides = enable;
            }

            public override string ToString()
            {
                var data = new List<(string key, string value)>()
                { 
                      ( nameof(pse_aac_sequence_classification_data         ),  pse_aac_sequence_classification_data           .ToString()     ) ,
                      ( nameof(sable_classification_data                    ),  sable_classification_data                      .ToString()     ) ,
                      ( nameof(mpsa_classification_data_subsequence         ),  mpsa_classification_data_subsequence           .ToString()     ) ,
                      ( nameof(blast_pssm_subsequence_classification_data   ),  blast_pssm_subsequence_classification_data     .ToString()     ) ,
                      ( nameof(aa_index_classification_data                 ),  aa_index_classification_data                   .ToString()     ) ,
                      ( nameof(sequence_geometry_classification_data        ),  sequence_geometry_classification_data          .ToString()     ) ,
                      ( nameof(intrinsically_unordered_data                 ),  intrinsically_unordered_data                   .ToString()     ) ,
                      ( nameof(dna_binding_prediction_data                  ),  dna_binding_prediction_data                    .ToString()     ) ,
                      ( nameof(r_peptides                                   ),  r_peptides                                     .ToString()     ) 
                };

                var ret = string.Join("\r\n", data.Select(a => a.key + " = " + a.value).ToList());

                return ret;
            }
        }

        public class feature_types_3d
        {
            public bool pse_ssc_dssp_classification_data = false;
            //public bool dssp_dist_classification_data = false;
            public bool foldx_classification_data = false;
            public bool ring_classification_data = false;
            public bool sasa_classification_data = false;
            public bool tortuosity_classification_data = false;
            public bool intramolecular_classification_data = false;
            public bool aa_aa_distances = false;

            public feature_types_3d()
            {

            }
            public feature_types_3d(bool enable)
            {
                pse_ssc_dssp_classification_data = enable;
                //dssp_dist_classification_data = enable;
                foldx_classification_data = enable;
                ring_classification_data = enable;
                sasa_classification_data = enable;
                tortuosity_classification_data = enable;
                intramolecular_classification_data = enable;
                aa_aa_distances = enable;
            }

            public override string ToString()
            {
                var data = new List<(string key, string value)>()
                {
                    ( nameof(pse_ssc_dssp_classification_data         ),  pse_ssc_dssp_classification_data           .ToString()     ) ,
                    ( nameof(foldx_classification_data         ),  foldx_classification_data           .ToString()     ) ,
                    ( nameof(ring_classification_data         ),  ring_classification_data           .ToString()     ) ,
                    ( nameof(sasa_classification_data         ),  sasa_classification_data           .ToString()     ) ,
                    ( nameof(tortuosity_classification_data         ),  tortuosity_classification_data           .ToString()     ) ,
                    ( nameof(intramolecular_classification_data         ),  intramolecular_classification_data           .ToString()     ) ,
                    ( nameof(aa_aa_distances         ),  aa_aa_distances           .ToString()     ) ,
                };

                var ret = string.Join("\r\n", data.Select(a => a.key + " = " + a.value).ToList());

                return ret;
            }
        }

        public static bool check_headers(List<feature_info> feats)
        {
            var header_list_str_dupe_check = feats.Select((a, i) => $"{a.alphabet},{a.dimension},{a.category},{a.source},{a.@group},{a.member},{a.perspective}").ToList();
            var header_list_str_dupe_check_distinct = header_list_str_dupe_check.Distinct().ToList();

            return header_list_str_dupe_check_distinct.Count == header_list_str_dupe_check.Count;
        }

        public static List<feature_info> calculate_classification_data_1d(subsequence_classification_data scd, List<Atom> subsequence_master_atoms, protein_data_sources source, int max_features, feature_types_1d feature_types_1d)
        {
#if DEBUG
            if (Program.verbose_debug) Program.WriteLine($"{nameof(calculate_classification_data_1d)}(subsequence_classification_data scd, List<Atom> subsequence_master_atoms, protein_data_sources source, int max_features, feature_types_1d feature_types_1d);");
#endif

            //var features_1d = new List<feature_info>();

            var features = new List<feature_info>();

            var check_num_features_consistency = true;

            var tasks = new List<Task<List<feature_info>>>();

            if (feature_types_1d != null)
            {
                if (feature_types_1d.pse_aac_sequence_classification_data)
                {
                    var task = Task.Run(() =>
                    {
                        var aa_seq_pse_aac_options = new pse_aac_options()
                        {
                            oaac = true,
                            oaac_binary = true,
                            motifs = true,
                            motifs_binary = true,
                            dipeptides = true,
                            dipeptides_binary = true,
                            //saac = true,
                            //saac_binary = true,
                            average_seq_position = true,
                            average_dipeptide_distance = true,
                        };

                        var pse_aac_sequence_classification_data = calculate_aa_or_ss_sequence_classification_data(source, "aa", "aa", scd.aa_subsequence, feature_calcs.seq_type.amino_acid_sequence, aa_seq_pse_aac_options, max_features);

                        
                        if (!check_headers(pse_aac_sequence_classification_data)) throw new Exception("duplicate headers");

                        pse_aac_sequence_classification_data = pse_aac_sequence_classification_data.GroupBy(a => (a.source, a.alphabet, a.category, a.dimension, a.@group)).Where(a => a.Count() <= max_features).SelectMany(a => a).ToList();

                        if (check_num_features_consistency)
                        {
                            var total_pse_aac_sequence_classification_data = pse_aac_sequence_classification_data?.Count ?? -1;
                            if (total_pse_aac_sequence_classification_data > -1 && subsequence_classification_data.total_pse_aac_sequence_classification_data > -1 && subsequence_classification_data.total_pse_aac_sequence_classification_data != total_pse_aac_sequence_classification_data) throw new Exception();
                            if (total_pse_aac_sequence_classification_data > -1) subsequence_classification_data.total_pse_aac_sequence_classification_data = total_pse_aac_sequence_classification_data;
                        }

                        return pse_aac_sequence_classification_data;
                    });
                    tasks.Add(task);
                }
                //

                if (feature_types_1d.sable_classification_data)
                {
                    var task = Task.Run(() =>
                    {
                        var sable_sequence_classification_data = calculate_sable_sequence_classification_data(scd, subsequence_master_atoms, source, max_features);
                        //features_1d.AddRange(sable_sequence_classification_data);

                        if (!check_headers(sable_sequence_classification_data)) throw new Exception("duplicate headers");

                        sable_sequence_classification_data = sable_sequence_classification_data.GroupBy(a => (a.source, a.alphabet, a.category, a.dimension, a.@group)).Where(a => a.Count() <= max_features).SelectMany(a => a).ToList();

                        if (check_num_features_consistency)
                        {
                            var total_sable_sequence_classification_data = sable_sequence_classification_data?.Count ?? -1;
                            if (total_sable_sequence_classification_data > -1 && subsequence_classification_data.total_sable_sequence_classification_data > -1 && subsequence_classification_data.total_sable_sequence_classification_data != total_sable_sequence_classification_data) throw new Exception();
                            if (total_sable_sequence_classification_data > -1) subsequence_classification_data.total_sable_sequence_classification_data = total_sable_sequence_classification_data;
                        }

                        return sable_sequence_classification_data;
                    }); tasks.Add(task);
                }

                if (feature_types_1d.mpsa_classification_data_subsequence)
                {
                    var task = Task.Run(() =>
                    {
                        var mpsa_classification_data = calculate_mpsa_classification_data(subsequence_master_atoms, source, max_features);

                        if (!check_headers(mpsa_classification_data)) throw new Exception("duplicate headers");

                        mpsa_classification_data = mpsa_classification_data.GroupBy(a => (a.source, a.alphabet, a.category, a.dimension, a.@group)).Where(a => a.Count() <= max_features).SelectMany(a => a).ToList();

                        

                        if (check_num_features_consistency)
                        {
                            var total_mpsa_classification_data = mpsa_classification_data?.Count ?? -1;
                            if (total_mpsa_classification_data > -1 && subsequence_classification_data.total_mpsa_classification_data > -1 && subsequence_classification_data.total_mpsa_classification_data != total_mpsa_classification_data) throw new Exception();
                            if (total_mpsa_classification_data > -1) subsequence_classification_data.total_mpsa_classification_data = total_mpsa_classification_data;

                        }

                        return mpsa_classification_data;
                    }); tasks.Add(task);
                }

                if (feature_types_1d.blast_pssm_subsequence_classification_data)
                {
                    var task = Task.Run(() =>
                    {
                        var pssm_options = new blast_pssm_options();

                        var blast_pssm_subsequence_classification_data = calculate_blast_pssm_classification_data(pssm_options, scd, subsequence_master_atoms, source, max_features);


                        if (!check_headers(blast_pssm_subsequence_classification_data)) throw new Exception("duplicate headers");

                        blast_pssm_subsequence_classification_data = blast_pssm_subsequence_classification_data.GroupBy(a => (a.source, a.alphabet, a.category, a.dimension, a.@group)).Where(a => a.Count() <= max_features).SelectMany(a => a).ToList();

                        if (check_num_features_consistency)
                        {
                            var total_blast_pssm_subsequence_classification_data = blast_pssm_subsequence_classification_data?.Count ?? -1;
                            if (total_blast_pssm_subsequence_classification_data > -1 && subsequence_classification_data.total_blast_pssm_subsequence_classification_data > -1 && subsequence_classification_data.total_blast_pssm_subsequence_classification_data != total_blast_pssm_subsequence_classification_data) throw new Exception();
                            if (total_blast_pssm_subsequence_classification_data > -1) subsequence_classification_data.total_blast_pssm_subsequence_classification_data = total_blast_pssm_subsequence_classification_data;
                        }

                        return blast_pssm_subsequence_classification_data;
                    }); tasks.Add(task);
                }

                if (feature_types_1d.aa_index_classification_data)
                {
                    var task = Task.Run(() =>
                    {
                        var aa_index_classification_data = calculate_aa_index_classification_data(subsequence_master_atoms, source, max_features);
                        //features_1d.AddRange(aa_index_classification_data);


                        if (!check_headers(aa_index_classification_data)) throw new Exception("duplicate headers");

                        aa_index_classification_data = aa_index_classification_data.GroupBy(a => (a.source, a.alphabet, a.category, a.dimension, a.@group)).Where(a => a.Count() <= max_features).SelectMany(a => a).ToList();

                        if (check_num_features_consistency)
                        {
                            var total_aa_index_classification_data = aa_index_classification_data?.Count ?? -1;
                            if (total_aa_index_classification_data > -1 && subsequence_classification_data.total_aa_index_classification_data > -1 && subsequence_classification_data.total_aa_index_classification_data != total_aa_index_classification_data) throw new Exception();
                            if (total_aa_index_classification_data > -1) subsequence_classification_data.total_aa_index_classification_data = total_aa_index_classification_data;
                        }

                        return aa_index_classification_data;
                    }); tasks.Add(task);
                }

                if (feature_types_1d.sequence_geometry_classification_data)
                {
                    var task = Task.Run(() =>
                    {
                        var sequence_geometry_classification_data = calculate_sequence_geometry_classification_data(scd, subsequence_master_atoms, source, max_features);
                        //features_1d.AddRange(sequence_geometry_classification_data);


                        if (!check_headers(sequence_geometry_classification_data)) throw new Exception("duplicate headers");

                        sequence_geometry_classification_data = sequence_geometry_classification_data.GroupBy(a => (a.source, a.alphabet, a.category, a.dimension, a.@group)).Where(a => a.Count() <= max_features).SelectMany(a => a).ToList();


                        if (check_num_features_consistency)
                        {
                            var total_sequence_geometry_classification_data = sequence_geometry_classification_data?.Count ?? -1;
                            if (total_sequence_geometry_classification_data > -1 && subsequence_classification_data.total_sequence_geometry_classification_data > -1 && subsequence_classification_data.total_sequence_geometry_classification_data != total_sequence_geometry_classification_data) throw new Exception();
                            if (total_sequence_geometry_classification_data > -1) subsequence_classification_data.total_sequence_geometry_classification_data = total_sequence_geometry_classification_data;
                        }

                        return sequence_geometry_classification_data;
                    }); tasks.Add(task);
                }

                if (feature_types_1d.intrinsically_unordered_data)
                {
                    var task = Task.Run(() =>
                    {
                        var intrinsically_unordered_data = calculate_intrinsically_unordered_data(subsequence_master_atoms, source, max_features);
                        //features_1d.AddRange(intrinsically_unordered_data);


                        if (!check_headers(intrinsically_unordered_data)) throw new Exception("duplicate headers");

                        intrinsically_unordered_data = intrinsically_unordered_data.GroupBy(a => (a.source, a.alphabet, a.category, a.dimension, a.@group)).Where(a => a.Count() <= max_features).SelectMany(a => a).ToList();

                        if (check_num_features_consistency)
                        {
                            var total_intrinsically_unordered_data = intrinsically_unordered_data?.Count ?? -1;
                            if (total_intrinsically_unordered_data > -1 && subsequence_classification_data.total_intrinsically_unordered_data > -1 && subsequence_classification_data.total_intrinsically_unordered_data != total_intrinsically_unordered_data) throw new Exception();
                            if (total_intrinsically_unordered_data > -1) subsequence_classification_data.total_intrinsically_unordered_data = total_intrinsically_unordered_data;
                        }

                        return intrinsically_unordered_data;
                    }); tasks.Add(task);
                }


                if (feature_types_1d.dna_binding_prediction_data)
                {
                    var task = Task.Run(() =>
                    {
                        var dna_binding_prediction_data = calculate_chain_dna_binding_prediction_data(subsequence_master_atoms, source, max_features);
                        //features_1d.AddRange(dna_binding_prediction_data);


                        if (!check_headers(dna_binding_prediction_data)) throw new Exception("duplicate headers");

                        dna_binding_prediction_data = dna_binding_prediction_data.GroupBy(a => (a.source, a.alphabet, a.category, a.dimension, a.@group)).Where(a => a.Count() <= max_features).SelectMany(a => a).ToList();

                        if (dna_binding_prediction_data != null && dna_binding_prediction_data.Count > 0)
                        {
                            if (check_num_features_consistency)
                            {
                                var total_dna_binding_prediction_data = dna_binding_prediction_data?.Count ?? -1;
                                if (total_dna_binding_prediction_data > -1 && subsequence_classification_data.total_dna_binding_prediction_data > -1 && subsequence_classification_data.total_dna_binding_prediction_data != total_dna_binding_prediction_data) throw new Exception();
                                if (total_dna_binding_prediction_data > -1) subsequence_classification_data.total_dna_binding_prediction_data = total_dna_binding_prediction_data;
                            }
                        }

                        return dna_binding_prediction_data;
                    }); tasks.Add(task);
                }


                if (feature_types_1d.r_peptides)
                {
                    var task = Task.Run(() =>
                    {
                        var seq = scd.aa_subsequence;

                        //var pep = new r_peptides();

                        var r_peptides_data = peptides_data(seq); //r_peptides.get_values(seq);

                        //try { pep.Dispose(); } catch (Exception) { } finally { }

                        r_peptides_data.ForEach(a =>
                        {
                            a.source = source.ToString();
                            a.alphabet = "Overall";
                        });

                        
                        if (!check_headers(r_peptides_data)) throw new Exception("duplicate headers");

                        //features_1d.AddRange(r_peptides_data);

                        r_peptides_data = r_peptides_data.GroupBy(a => (a.source, a.alphabet, a.category, a.dimension, a.@group)).Where(a => a.Count() <= max_features).SelectMany(a => a).ToList();

                        if (check_num_features_consistency)
                        {
                            var total_r_peptides_prediction_data = r_peptides_data?.Count ?? -1;
                            if (total_r_peptides_prediction_data > -1 && subsequence_classification_data.total_r_peptides_prediction_data > -1 && subsequence_classification_data.total_r_peptides_prediction_data != total_r_peptides_prediction_data) throw new Exception();
                            if (total_r_peptides_prediction_data > -1) subsequence_classification_data.total_r_peptides_prediction_data = total_r_peptides_prediction_data;
                        }

                        return r_peptides_data;

                    });
                    tasks.Add(task);
                }
            }

            Task.WaitAll(tasks.ToArray<Task>());

            tasks.ForEach(a => features.AddRange(a.Result));







            return features;
        }

        
        private static int peptides_call_count = 0;
        private static object peptides_call_count_lock = new object();

        public static List<feature_info> peptides_data(string sequence)
        {
            var call_count = 0;
            lock (peptides_call_count_lock)
            {
                call_count = peptides_call_count++;
            }
#if DEBUG
            var exe = @"C:\Users\aaron\Desktop\Projects\Coil_DHC_Dataset_Maker\peptides_server\bin\x64\Debug\peptides_server.exe";
#else
            var exe = @"C:\Users\aaron\Desktop\Projects\Coil_DHC_Dataset_Maker\peptides_server\bin\x64\Release\peptides_server.exe";
#endif

            //var psi = new ProcessStartInfo()
            //{
            //    Arguments = sequence,
            //    CreateNoWindow = false,
            //    FileName = exe,
            //    WorkingDirectory = Path.GetDirectoryName(exe),
            //    RedirectStandardInput = true,
            //    RedirectStandardError = true,
            //    UseShellExecute = false,
            //    WindowStyle = ProcessWindowStyle.Hidden,
            //};

            var start = new ProcessStartInfo
            {
                FileName = exe,
                Arguments = $"{call_count} {sequence}",
                UseShellExecute = false,
                CreateNoWindow = false,
                RedirectStandardOutput = true,
                RedirectStandardError = true,
                //WindowStyle = ProcessWindowStyle.Hidden,

            };


            var result = new List<feature_info>();

            using (var process = Process.Start(start))
            {
                process.PriorityBoostEnabled = true;
                process.PriorityClass = ProcessPriorityClass.High;


                using (var reader = process.StandardOutput)
                {
                    var data = reader.ReadToEnd();
                    var stderr = process.StandardError.ReadToEnd();

                    data = data.Substring(data.IndexOf('\n') + 1);



                    process.WaitForExit();

                    //Console.WriteLine("Data: " + data);

                    var r = feature_info_container.deserialise(data);

                    result = r.feautre_info_list;
                }
            }

            return result;
        }

        public static List<feature_info> calculate_classification_data_3d(subsequence_classification_data scd, List<Atom> subsequence_master_atoms, protein_data_sources source, int max_features, feature_types_3d feature_types_3d)
        {
#if DEBUG
            if (Program.verbose_debug) Program.WriteLine($"{nameof(calculate_classification_data_3d)}(subsequence_classification_data scd, List<Atom> subsequence_master_atoms, protein_data_sources source, int max_features, feature_types_3d feature_types_3d);");
#endif

            var tasks = new List<Task<List<feature_info>>>();

            var features = new List<feature_info>();

            var check_num_features_consistency = true;

            if (feature_types_3d != null)
            {
                //if (feature_types_3d.dssp_dist_classification_data)
                //{
                //    var make_dssp_feature = true;
                //    var make_stride_feature = false;

                //    var protein_dssp_dist_classification_data = calculate_dssp_and_stride_protein_classification_data(make_dssp_feature, make_stride_feature, source, scd, max_features);

                //    features_3d.AddRange(protein_dssp_dist_classification_data);

                //    if (check_num_features_consistency)
                //    {
                //        var total_protein_dssp_dist_classification_data = protein_dssp_dist_classification_data?.Count ?? -1;
                //        if (total_protein_dssp_dist_classification_data > -1 && subsequence_classification_data.total_protein_dssp_dist_classification_data > -1 && subsequence_classification_data.total_protein_dssp_dist_classification_data != total_protein_dssp_dist_classification_data) throw new Exception();
                //        if (total_protein_dssp_dist_classification_data > -1) subsequence_classification_data.total_protein_dssp_dist_classification_data = total_protein_dssp_dist_classification_data;
                //    }
                //}

                if (feature_types_3d.pse_ssc_dssp_classification_data)
                {
                    var task = Task.Run(() =>
                    {
                        var aa_seq_pse_aac_options = new pse_aac_options()
                        {
                            oaac = true,
                            oaac_binary = true,
                            motifs = true,
                            motifs_binary = true,
                            dipeptides = true,
                            dipeptides_binary = true,
                            //saac = true,
                            //saac_binary = true,
                            average_seq_position = true,
                            average_dipeptide_distance = true,
                        };

                        var pse_ssc_dssp_classification_data = calculate_aa_or_ss_sequence_classification_data(source, "dssp_monomer", "dssp_monomer", scd.dssp_monomer_subsequence, feature_calcs.seq_type.secondary_structure_sequence, aa_seq_pse_aac_options, max_features);


                        if (!check_headers(pse_ssc_dssp_classification_data)) throw new Exception("duplicate headers");

                        pse_ssc_dssp_classification_data = pse_ssc_dssp_classification_data.GroupBy(a => (a.source, a.alphabet, a.category, a.dimension, a.@group)).Where(a => a.Count() <= max_features).SelectMany(a => a).ToList();

                        if (check_num_features_consistency)
                        {
                            var total_pse_ssc_dssp_classification_data = pse_ssc_dssp_classification_data?.Count ?? -1;
                            if (total_pse_ssc_dssp_classification_data > -1 && subsequence_classification_data.total_pse_ssc_dssp_classification_data > -1 && subsequence_classification_data.total_pse_ssc_dssp_classification_data != total_pse_ssc_dssp_classification_data) throw new Exception();
                            if (total_pse_ssc_dssp_classification_data > -1) subsequence_classification_data.total_pse_ssc_dssp_classification_data = total_pse_ssc_dssp_classification_data;
                        }

                        return pse_ssc_dssp_classification_data;
                    });
                    tasks.Add(task);
                }

                if (feature_types_3d.foldx_classification_data)
                {
                    var task = Task.Run(() =>
                    {
                        var foldx_classification_data = calculate_foldx_classification_data(scd, subsequence_master_atoms, source, max_features);


                        if (!check_headers(foldx_classification_data)) throw new Exception("duplicate headers");

                        foldx_classification_data = foldx_classification_data.GroupBy(a => (a.source, a.alphabet, a.category, a.dimension, a.@group)).Where(a => a.Count() <= max_features).SelectMany(a => a).ToList();

                        

                        if (check_num_features_consistency)
                        {
                            if (source == protein_data_sources.subsequence_3d)
                            {
                                var total_foldx_classification_subsequence_3d_data = foldx_classification_data?.Count ?? -1;
                                if (total_foldx_classification_subsequence_3d_data > -1 && subsequence_classification_data.total_foldx_classification_subsequence_3d_data > -1 && subsequence_classification_data.total_foldx_classification_subsequence_3d_data != total_foldx_classification_subsequence_3d_data) throw new Exception();
                                if (total_foldx_classification_subsequence_3d_data > -1) subsequence_classification_data.total_foldx_classification_subsequence_3d_data = total_foldx_classification_subsequence_3d_data;
                            }
                            else if (source == protein_data_sources.neighbourhood_3d)
                            {
                                var total_foldx_classification_neighbourhood_3d_data = foldx_classification_data?.Count ?? -1;
                                if (total_foldx_classification_neighbourhood_3d_data > -1 && subsequence_classification_data.total_foldx_classification_neighbourhood_3d_data > -1 && subsequence_classification_data.total_foldx_classification_neighbourhood_3d_data != total_foldx_classification_neighbourhood_3d_data) throw new Exception();
                                if (total_foldx_classification_neighbourhood_3d_data > -1) subsequence_classification_data.total_foldx_classification_neighbourhood_3d_data = total_foldx_classification_neighbourhood_3d_data;
                            }
                            else if (source == protein_data_sources.protein_3d)
                            {
                                var total_foldx_classification_protein_3d_data = foldx_classification_data?.Count ?? -1;
                                if (total_foldx_classification_protein_3d_data > -1 && subsequence_classification_data.total_foldx_classification_protein_3d_data > -1 && subsequence_classification_data.total_foldx_classification_protein_3d_data != total_foldx_classification_protein_3d_data) throw new Exception();
                                if (total_foldx_classification_protein_3d_data > -1) subsequence_classification_data.total_foldx_classification_protein_3d_data = total_foldx_classification_protein_3d_data;
                            }

                        }

                        return foldx_classification_data;
                    });
                    tasks.Add(task);
                }

                if (feature_types_3d.ring_classification_data)
                {
                    var task = Task.Run(() =>
                    {
                        var ring_classification_data = calculate_ring_classification_data(scd, subsequence_master_atoms, source, max_features);


                        if (!check_headers(ring_classification_data)) throw new Exception("duplicate headers");

                        ring_classification_data = ring_classification_data.GroupBy(a => (a.source, a.alphabet, a.category, a.dimension, a.@group)).Where(a => a.Count() <= max_features).SelectMany(a => a).ToList();

                        if (check_num_features_consistency)
                        {
                            var total_ring_classification_data = ring_classification_data?.Count ?? -1;
                            if (total_ring_classification_data > -1 && subsequence_classification_data.total_ring_classification_data > -1 && subsequence_classification_data.total_ring_classification_data != total_ring_classification_data) throw new Exception();
                            if (total_ring_classification_data > -1) subsequence_classification_data.total_ring_classification_data = total_ring_classification_data;
                        }

                        return ring_classification_data;
                    });
                    tasks.Add(task);
                }

                if (feature_types_3d.sasa_classification_data)
                {
                    var task = Task.Run(() =>
                    {
                        var sasa_classification_data = calculate_sasa_classification_data(subsequence_master_atoms, source, max_features);


                        if (!check_headers(sasa_classification_data)) throw new Exception("duplicate headers");

                        sasa_classification_data = sasa_classification_data.GroupBy(a => (a.source, a.alphabet, a.category, a.dimension, a.@group)).Where(a => a.Count() <= max_features).SelectMany(a => a).ToList();

                        if (check_num_features_consistency)
                        {
                            var total_sasa_classification_data = sasa_classification_data?.Count ?? -1;
                            if (total_sasa_classification_data > -1 && subsequence_classification_data.total_sasa_classification_data > -1 && subsequence_classification_data.total_sasa_classification_data != total_sasa_classification_data) throw new Exception();
                            if (total_sasa_classification_data > -1) subsequence_classification_data.total_sasa_classification_data = total_sasa_classification_data;

                        }

                        return sasa_classification_data;
                    });
                    tasks.Add(task);

                }

                if (feature_types_3d.tortuosity_classification_data)
                {
                    var task = Task.Run(() =>
                    {
                        var tortuosity_classification_data = calculate_tortuosity_classification_data(scd, subsequence_master_atoms, source, max_features);


                        if (!check_headers(tortuosity_classification_data)) throw new Exception("duplicate headers");

                        tortuosity_classification_data = tortuosity_classification_data.GroupBy(a => (a.source, a.alphabet, a.category, a.dimension, a.@group)).Where(a => a.Count() <= max_features).SelectMany(a => a).ToList();

                        if (check_num_features_consistency)
                        {
                            var total_tortuosity_classification_data = tortuosity_classification_data?.Count ?? -1;
                            if (total_tortuosity_classification_data > -1 && subsequence_classification_data.total_tortuosity_classification_data > -1 && subsequence_classification_data.total_tortuosity_classification_data != total_tortuosity_classification_data) throw new Exception();
                            if (total_tortuosity_classification_data > -1) subsequence_classification_data.total_tortuosity_classification_data = total_tortuosity_classification_data;

                        }

                        return tortuosity_classification_data;
                    });
                    tasks.Add(task);
                }

                if (feature_types_3d.intramolecular_classification_data)
                {
                    var task = Task.Run(() =>
                    {
                        var intramolecular_classification_data = calculate_intramolecular_classification_data(subsequence_master_atoms, source, max_features);

                        if (!check_headers(intramolecular_classification_data)) throw new Exception("duplicate headers");

                        intramolecular_classification_data = intramolecular_classification_data.GroupBy(a => (a.source, a.alphabet, a.category, a.dimension, a.@group)).Where(a => a.Count() <= max_features).SelectMany(a => a).ToList();

                        if (check_num_features_consistency)
                        {
                            var total_intramolecular_classification_data = intramolecular_classification_data?.Count ?? -1;
                            if (total_intramolecular_classification_data > -1 && subsequence_classification_data.total_intramolecular_classification_data > -1 && subsequence_classification_data.total_intramolecular_classification_data != total_intramolecular_classification_data) throw new Exception();
                            if (total_intramolecular_classification_data > -1) subsequence_classification_data.total_intramolecular_classification_data = total_intramolecular_classification_data;
                        }

                        return intramolecular_classification_data;
                    });
                    tasks.Add(task);
                }

                if (feature_types_3d.aa_aa_distances)
                {
                    var task = Task.Run(() =>
                    {
                        var aa_aa_distances_classification_data = calculate_aa_aa_distances_classification_data(subsequence_master_atoms, source, max_features);

                        if (!check_headers(aa_aa_distances_classification_data)) throw new Exception("duplicate headers");

                        aa_aa_distances_classification_data = aa_aa_distances_classification_data.GroupBy(a => (a.source, a.alphabet, a.category, a.dimension, a.@group)).Where(a => a.Count() <= max_features).SelectMany(a => a).ToList();

                        if (check_num_features_consistency)
                        {
                            var total_aa_aa_distances_classification_data = aa_aa_distances_classification_data?.Count ?? -1;
                            if (total_aa_aa_distances_classification_data > -1 && subsequence_classification_data.total_aa_aa_distances_classification_data > -1 && subsequence_classification_data.total_aa_aa_distances_classification_data != total_aa_aa_distances_classification_data) throw new Exception();
                            if (total_aa_aa_distances_classification_data > -1) subsequence_classification_data.total_aa_aa_distances_classification_data = total_aa_aa_distances_classification_data;
                        }

                        return aa_aa_distances_classification_data;
                    });
                    tasks.Add(task);
                }
            }

            Task.WaitAll(tasks.ToArray<Task>());

            tasks.ForEach(a => features.AddRange(a.Result));

            return features;
        }


        //public static List<feature_info> calculate_dssp_and_stride_subsequence_classification_data_template = null;

        //public static List<feature_info> calculate_dssp_and_stride_subsequence_classification_data(bool make_dssp_feature, bool make_stride_feature, protein_data_sources source, subsequence_classification_data subsequence_classification_data, int max_features)
        //{
        //    if (subsequence_classification_data.subsequence_master_atoms == null || subsequence_classification_data.subsequence_master_atoms.Count == 0)
        //    {
        //        var template = calculate_dssp_and_stride_subsequence_classification_data_template.Select(a => new feature_info(a)
        //        {
        //            source = source,
        //            feature_value = 0
        //        }).ToList();

        //        return template;
        //    }

        //    //if (make_dssp_feature || make_stride_feature)
        //    //{
        //    //    var x = encode_dssp_stride(make_dssp_feature, make_stride_feature, source, subsequence_classification_data, max_features);
        //    //    row_features.AddRange(x);
        //    //}

        //    // sec struct (should only be CCCCCCCCCC for coils; EEEEEEEEE for dimophics/dhc) - useful to test their neighbourhoods though?
        //    // actually not, since we observed some variation in Coil/Strand between chain A and B
        //    // updated to monomeric structures, not multimeric
        //    var row_features = new List<feature_info>();


        //    for (var count_or_dist = 0; count_or_dist <= 1; count_or_dist++)
        //    {
        //        for (var normal_or_sqrt = 0; normal_or_sqrt <= 1; normal_or_sqrt++)
        //        {
        //            var as_sqrt = normal_or_sqrt != 0;
        //            var as_dist = count_or_dist != 0;

        //            if (!as_dist) continue;

        //            var dist_name = $"{(as_dist ? "dist" : "count")}_{(as_sqrt ? "sqrt" : "normal")}";


        //            foreach (feature_calcs.ss_types ss_type in Enum.GetValues(typeof(feature_calcs.ss_types)))
        //            {
        //                var ss_type_name = Enum.GetName(typeof(feature_calcs.ss_types), ss_type)?.ToLowerInvariant();
        //                string ss_seq = "";

        //                if (ss_type == feature_calcs.ss_types.DSSP)
        //                {
        //                    ss_seq = subsequence_classification_data.dssp_monomer_subsequence;
        //                }
        //                else if (ss_type == feature_calcs.ss_types.DSSP3)
        //                {
        //                    ss_seq = subsequence_classification_data.dssp_monomer_subsequence;
        //                }
        //                else if (ss_type == feature_calcs.ss_types.STRIDE)
        //                {
        //                    ss_seq = subsequence_classification_data.stride_monomer_subsequence;
        //                }
        //                else if (ss_type == feature_calcs.ss_types.STRIDE3)
        //                {
        //                    ss_seq = subsequence_classification_data.stride_monomer_subsequence;
        //                }

        //                if (make_dssp_feature && (ss_type == feature_calcs.ss_types.DSSP || ss_type == feature_calcs.ss_types.DSSP3))
        //                {
        //                    var ss_distribution = feature_calcs.ss_distribution(ss_seq, as_sqrt, as_dist, ss_type);
        //                    var xf = ss_distribution.Select(a => { return new feature_info() { alphabet = "Overall", dimension = 3, category = "dssp", source = source, @group = $"{nameof(ss_distribution)}_{ss_type_name}_{dist_name}", member = a.ss_type.ToString(), perspective = $@"default", feature_value = a.value }; }).ToList();
        //                    if (xf.Count <= max_features) row_features.AddRange(xf);
        //                }

        //                if (make_stride_feature && (ss_type == feature_calcs.ss_types.STRIDE || ss_type == feature_calcs.ss_types.STRIDE3))
        //                {
        //                    var ss_distribution = feature_calcs.ss_distribution(ss_seq, as_sqrt, as_dist, ss_type);
        //                    var xf = ss_distribution.Select(a => new feature_info() { alphabet = "Overall", dimension = 3, source = source, category = "stride", group = $"{nameof(ss_distribution)}_{ss_type_name}_{dist_name}", member = a.ss_type.ToString(), perspective = $@"default", feature_value = a.value }).ToList();
        //                    if (xf.Count <= max_features) row_features.AddRange(xf);

        //                }
        //            }
        //        }
        //    }

        //    if (calculate_dssp_and_stride_subsequence_classification_data_template == null)
        //    {
        //        var template = row_features.Select(a => new feature_info(a) { source = "", feature_value = 0 }).ToList();
        //        calculate_dssp_and_stride_subsequence_classification_data_template = template;
        //    }

        //    return row_features;

        //}

//        public static List<feature_info> calculate_dssp_and_stride_protein_classification_data_template = null;

//        public static List<feature_info> calculate_dssp_and_stride_protein_classification_data(bool make_dssp_feature, bool make_stride_feature, protein_data_sources source, subsequence_classification_data subsequence_classification_data, int max_features)
//        {
//#if DEBUG
//            if (Program.verbose_debug) Program.WriteLine($"{nameof(calculate_dssp_and_stride_protein_classification_data)}(bool make_dssp_feature, bool make_stride_feature, protein_data_sources source, subsequence_classification_data subsequence_classification_data, int max_features);");
//#endif

//            if (subsequence_classification_data.subsequence_master_atoms == null || subsequence_classification_data.subsequence_master_atoms.Count == 0)
//            {
//                if (calculate_dssp_and_stride_protein_classification_data_template == null) throw new Exception();

//                var template = calculate_dssp_and_stride_protein_classification_data_template.Select(a => new feature_info(a)
//                {
//                    source = source.ToString(),
//                    feature_value = 0
//                }).ToList();

//                return template;
//            }
//            var row_features = new List<feature_info>();

//            //var chain_atoms = subsequence_classification_data.subsequence_master_atoms.Where(a => a != null && a.chain_atoms != null && a.chain_atoms.Count > 0).FirstOrDefault()?.chain_atoms;
//            //var chain_master_atoms = Atom.select_amino_acid_master_atoms(null, chain_atoms);

//            //var stride_seq = string.Join("", chain_master_atoms.Select(a => a.monomer_stride).ToList());
//            //var dssp_seq = string.Join("", chain_master_atoms.Select(a => a.monomer_dssp).ToList());

//            var dssp_seq = subsequence_classification_data.dssp_monomer_subsequence;
//            var stride_seq = subsequence_classification_data.stride_monomer_subsequence;

//            for (var count_or_dist = 0; count_or_dist <= 1; count_or_dist++)
//            {
//                for (var normal_or_sqrt = 0; normal_or_sqrt <= 1; normal_or_sqrt++)
//                {
//                    var as_sqrt = normal_or_sqrt != 0;
//                    var as_dist = count_or_dist != 0;

//                    if (!as_dist) continue;

//                    var dist_name = $"{(as_dist ? "dist" : "count")}_{(as_sqrt ? "sqrt" : "normal")}";


//                    foreach (feature_calcs.ss_types ss_type in Enum.GetValues(typeof(feature_calcs.ss_types)))
//                    {
//                        var ss_type_name = Enum.GetName(typeof(feature_calcs.ss_types), ss_type)?.ToLowerInvariant();
//                        string ss_seq = "";

//                        if (ss_type == feature_calcs.ss_types.DSSP)
//                        {
//                            ss_seq = dssp_seq;
//                        }
//                        else if (ss_type == feature_calcs.ss_types.DSSP3)
//                        {
//                            ss_seq = dssp_seq;
//                        }
//                        else if (ss_type == feature_calcs.ss_types.STRIDE)
//                        {
//                            ss_seq = stride_seq;
//                        }
//                        else if (ss_type == feature_calcs.ss_types.STRIDE3)
//                        {
//                            ss_seq = stride_seq;
//                        }

//                        if (make_dssp_feature && (ss_type == feature_calcs.ss_types.DSSP || ss_type == feature_calcs.ss_types.DSSP3))
//                        {
//                            var ss_distribution = feature_calcs.ss_distribution(ss_seq, as_sqrt, as_dist, ss_type);
//                            var ss_distribution_feats = ss_distribution.Select(a =>

//                                new feature_info()
//                                {
//                                    alphabet = "Overall",
//                                    dimension = 3,
//                                    category = "dssp_monomer",
//                                    source = source.ToString(),
//                                    @group = $"dssp_monomer_{nameof(ss_distribution)}_{ss_type_name}_{dist_name}",
//                                    member = a.ss_type.ToString(),
//                                    perspective = $@"default",
//                                    feature_value = a.value
//                                }
//                            ).ToList();

//                            if (ss_distribution_feats.Count <= max_features) row_features.AddRange(ss_distribution_feats);
//                        }

//                        if (make_stride_feature && (ss_type == feature_calcs.ss_types.STRIDE || ss_type == feature_calcs.ss_types.STRIDE3))
//                        {
//                            var ss_distribution = feature_calcs.ss_distribution(ss_seq, as_sqrt, as_dist, ss_type);
//                            var ss_distribution_feats = ss_distribution.Select(a =>
//                                new feature_info()
//                                {
//                                    alphabet = "Overall",
//                                    dimension = 3,
//                                    source = source.ToString(),
//                                    category = "stride_monomer",
//                                    group = $"stride_monomer_{nameof(ss_distribution)}_{ss_type_name}_{dist_name}",
//                                    member = a.ss_type.ToString(),
//                                    perspective = $@"default",
//                                    feature_value = a.value
//                                }
//                            ).ToList();
//                            if (ss_distribution_feats.Count <= max_features) row_features.AddRange(ss_distribution_feats);

//                        }
//                    }
//                }
//            }

//            if (calculate_dssp_and_stride_protein_classification_data_template == null)
//            {
//                var template = row_features.Select(a => new feature_info(a) { source = "", feature_value = 0 }).ToList();
//                calculate_dssp_and_stride_protein_classification_data_template = template;
//            }

//            return row_features;

//        }

        //private static bool IsBitSet(int b, int pos)
        //{
        //    return (b & (1 << pos)) != 0;
        //}

        public static List<feature_info> calculate_aa_or_ss_sequence_classification_data_template = null;

        public static List<feature_info> calculate_aa_or_ss_sequence_classification_data(protein_data_sources source, string category_prefix, string group_prefix, string sequence, feature_calcs.seq_type seq_type, pse_aac_options pse_aac_options, int max_features)
        {
#if DEBUG
            if (Program.verbose_debug) Program.WriteLine($"{nameof(calculate_aa_or_ss_sequence_classification_data)}(protein_data_sources source, string category_prefix, string group_prefix, string sequence, feature_calcs.seq_type seq_type, feature_calcs.pse_aac_options pse_aac_options, int max_features);");
#endif

            if (sequence == null || sequence.Length == 0)
            {
                if (calculate_aa_or_ss_sequence_classification_data_template == null) throw new Exception();

                var template = calculate_aa_or_ss_sequence_classification_data_template.Select(a => new feature_info(a)
                {
                    source = source.ToString(),
                    feature_value = 0
                }).ToList();

                return template;
            }

            var features = new List<feature_info>();

            for (var count_or_dist = 0; count_or_dist <= 1; count_or_dist++)
            {
                for (var normal_or_sqrt = 0; normal_or_sqrt <= 1; normal_or_sqrt++)
                {
                    var as_sqrt = normal_or_sqrt != 0;
                    var as_dist = count_or_dist != 0;

                    if (!as_dist) continue;// skip non-dist (i.e. total count) for now, as it doesn't seem to add much value (if any)
                    if (as_sqrt) continue; // skip sqrt for now, as it doesn't seem to add much value (if any)


                    var dist_name = $"{(as_dist ? "dist" : "count")}_{(as_sqrt ? "sqrt" : "normal")}";

                    var seq_split = feature_calcs.split_sequence(sequence, 3, false);

                    var seqs = new List<(string name, string sequence)>();

                    seqs.Add(("unsplit", sequence));
                    seqs.AddRange(seq_split.Select(a=>("split",a)).ToList());

                    for (var sq_index = 0; sq_index < seqs.Count; sq_index++)
                    {
                        var sq = seqs[sq_index];
                        var seq_feature_alphabets = feature_calcs.feature_pse_aac(sq.sequence, seq_type, pse_aac_options, as_sqrt, as_dist);

                        foreach (var seq_feature_alphabet in seq_feature_alphabets)
                        {
                            var ( alphabet_id, alphabet_name, motifs, motifs_binary,
                                //saac,
                                //saac_binary,
                                oaac, oaac_binary, average_seq_positions, dipeptides, dipeptides_binary, average_dipeptide_distance ) = seq_feature_alphabet;



                            if (pse_aac_options.dipeptides || pse_aac_options.dipeptides_binary)
                            {
                                var dipeptides_list = new List<(string name, feature_calcs.named_double[][] dipeptides)>();
                                if (pse_aac_options.dipeptides) dipeptides_list.Add((nameof(dipeptides), dipeptides));
                                if (pse_aac_options.dipeptides_binary) dipeptides_list.Add((nameof(dipeptides_binary), dipeptides_binary));

                                foreach (var d in dipeptides_list)
                                {
                                    // list [distance] [aa to aa]

                                    for (var distance = 0; distance < d.dipeptides.Length; distance++)
                                    {
                                        var f_list = new List<feature_info>();

                                        for (var aa_to_aa = 0; aa_to_aa < d.dipeptides[distance].Length; aa_to_aa++)
                                        {
                                            var item = d.dipeptides[distance][aa_to_aa];

                                            var f = new feature_info()
                                            {
                                                alphabet = alphabet_name,
                                                dimension = 1,
                                                category = $@"{category_prefix}_{d.name}",
                                                source = source.ToString(),
                                                @group = $@"{group_prefix}_{sq.name}_{d.name}_{(distance+1)}_{alphabet_name}_{dist_name}",
                                                member = $@"{sq_index}_{aa_to_aa}_{item.name}",
                                                perspective = $@"default",
                                                feature_value = item.value
                                            };

                                            f_list.Add(f);
                                        }

                                        if (f_list.Count <= max_features) features.AddRange(f_list);
                                    }

                                    /*

                                    var order_context_values_joined = new List<(string x, List<feature_calcs.named_double> y)>();
                                    var order_context_bits = d.dipeptides.Length;
                                    var order_context_max_joined_contexts = Convert.ToInt32(new string('1', d.dipeptides.Length), 2);
                                    for (var order_context_combin_id = 1; order_context_combin_id <= order_context_max_joined_contexts; order_context_combin_id++)
                                    {
                                        // join every possible combination of the contexts/distances

                                        var order_context_joined = new List<feature_calcs.named_double>();
                                        var order_context_distance_lengths = new List<int>();

                                        for (var order_context_bit_id = 0; order_context_bit_id < order_context_bits; order_context_bit_id++)
                                        {
                                            if (IsBitSet(order_context_combin_id, order_context_bit_id))
                                            {
                                                order_context_distance_lengths.Add(order_context_bit_id);
                                                order_context_joined.AddRange(d.dipeptides[order_context_bit_id]);
                                            }
                                        }

                                        order_context_values_joined.Add((string.Join("_", order_context_distance_lengths), order_context_joined));
                                    }
                                    

                                    for (var index = 0; index < order_context_values_joined.Count; index++)
                                    {
                                        // this adds all possible combinations of different lengths of order_context (001, 010, 100, 101, 110, 011, 111)
                                        var order_context_vj = order_context_values_joined[index];
                                        var x7 = order_context_vj.y.Select(x => new feature_info()
                                        {
                                            alphabet = alphabet_name,
                                            dimension = 1,
                                            category = $@"{category_prefix}_{d.name}",
                                            source = source.ToString(),
                                            @group = $@"{group_prefix}_{sq.name}_{d.name}_{order_context_vj.x}_{alphabet_name}_{dist_name}_{index}",
                                            member = $@"{sq_index}_{x.name}",
                                            perspective = $@"default",
                                            feature_value = x.value
                                        }).ToList();

                                        if (x7.Count <= max_features) features.AddRange(x7);
                                    }
                                    */


                                }
                            }
                            
                            //"Physicochemical,1,aa_dipeptides,subsequence_1d,aa_unsplit_dipeptides_0_1_2_Physicochemical_dist_normal_6,0_AVFPMILW_AVFPMILW,default"
                            //continue;

                            if (pse_aac_options.motifs || pse_aac_options.motifs_binary)
                            {
                                // motifs [motif length] [motif]

                                var motifs_list = new List<(string name, feature_calcs.named_double[][] motifs)>();
                                if (pse_aac_options.motifs) motifs_list.Add((nameof(motifs), motifs));
                                if (pse_aac_options.motifs_binary) motifs_list.Add((nameof(motifs_binary), motifs_binary));

                                foreach (var m in motifs_list)
                                {
                                    for (var motif_length = 0; motif_length < m.motifs.Length; motif_length++)
                                    {
                                        var f_list = new List<feature_info>();

                                        for (var aa_motif_index = 0; aa_motif_index < m.motifs[motif_length].Length; aa_motif_index++)
                                        {
                                            var item = m.motifs[motif_length][aa_motif_index];

                                            var f = new feature_info()
                                            {
                                                alphabet = alphabet_name,
                                                dimension = 1,
                                                category = $@"{category_prefix}_{m.name}",
                                                source = source.ToString(),
                                                @group = $@"{group_prefix}_{sq.name}_{m.name}_{(motif_length + 1)}_{alphabet_name}_{dist_name}",
                                                member = $@"{sq_index}_{aa_motif_index}_{item.name}",
                                                perspective = $@"default",
                                                feature_value = item.value
                                            };

                                            f_list.Add(f);
                                        }

                                        if (f_list.Count <= max_features) features.AddRange(f_list);
                                    }

                                    /*var motifs_values_joined = new List<(string x, List<feature_calcs.named_double> y)>();
                                    var motifs_bits = m.motifs.Length;
                                    var motifs_max_joined_contexts = Convert.ToInt32(new string('1', m.motifs.Length), 2);
                                    for (var motifs_combin_id = 1; motifs_combin_id <= motifs_max_joined_contexts; motifs_combin_id++)
                                    {
                                        // join every possible combination of the motifs/lengths

                                        var motifs_joined = new List<feature_calcs.named_double>();
                                        var motifs_distance_lengths = new List<int>();

                                        for (var motifs_bit_id = 0; motifs_bit_id < motifs_bits; motifs_bit_id++)
                                        {
                                            if (IsBitSet(motifs_combin_id, motifs_bit_id))
                                            {
                                                motifs_distance_lengths.Add(motifs_bit_id);
                                                motifs_joined.AddRange(m.motifs[motifs_bit_id]);
                                            }
                                        }

                                        motifs_values_joined.Add((string.Join("_", motifs_distance_lengths), motifs_joined));
                                    }

                                    for (var index = 0; index < motifs_values_joined.Count; index++)
                                    {
                                        // this adds all possible combinations of different lengths of motifs (001, 010, 100, 101, 110, 011, 111)
                                        var motifs_vj = motifs_values_joined[index];
                                        var x7 = motifs_vj.y.Select(x => new feature_info()
                                        {
                                            alphabet = alphabet_name,
                                            dimension = 1,
                                            category = $@"{category_prefix}_{m.name}",
                                            source = source.ToString(),
                                            @group = $@"{group_prefix}_{sq.name}_{m.name}_{motifs_vj.x}_{alphabet_name}_{dist_name}_{index}",
                                            member = $@"{sq_index}_{x.name}",
                                            perspective = $@"default",
                                            feature_value = x.value
                                        }).ToList();

                                        if (x7.Count <= max_features) features.AddRange(x7);
                                    }*/
                                }
                            }


                            //if (pse_aac_options.saac || pse_aac_options.saac_binary) //split_composition_values != null)// && split_composition_values.Length>0)
                            //{
                            //    var saac_list = new List<(string name, feature_calcs.named_double[] saac)>();
                            //    if (pse_aac_options.saac) saac_list.Add((nameof(saac), saac));
                            //    if (pse_aac_options.saac_binary) saac_list.Add((nameof(saac_binary), saac_binary));

                            //    foreach (var s in saac_list)
                            //    {
                            //        var x2 = s.saac.Select(x => new feature_info()
                            //        {
                            //            alphabet = alphabet_name,
                            //            dimension = 1,
                            //            category = $@"{category_prefix}_{s.name}",
                            //            source = source.ToString(),
                            //            group = $@"{group_prefix}_{s.name}_{alphabet_name}_{dist_name}",
                            //            member = $@"{sq_index}_{x.name}",
                            //            perspective = $@"default",
                            //            feature_value = x.value
                            //        }).ToList();

                            //        if (x2.Count <= max_features) features.AddRange(x2);
                            //    }
                            //}

                            if (pse_aac_options.oaac || pse_aac_options.oaac_binary) //composition_values != null)// && composition_values.Length>0)
                            {
                                var oaac_list = new List<(string name, feature_calcs.named_double[] oaac)>();
                                if (pse_aac_options.oaac) oaac_list.Add((nameof(oaac), oaac));
                                if (pse_aac_options.oaac_binary) oaac_list.Add((nameof(oaac_binary), oaac_binary));

                                foreach (var o in oaac_list)
                                {
                                    var x3 = o.oaac.Select(x => new feature_info()
                                    {
                                        alphabet = alphabet_name,
                                        dimension = 1,
                                        category = $@"{category_prefix}_{o.name}",
                                        source = source.ToString(),
                                        @group = $@"{group_prefix}_{sq.name}_{o.name}_{alphabet_name}_{dist_name}",
                                        member = $@"{sq_index}_{x.name}",
                                        perspective = $@"default",
                                        feature_value = x.value
                                    }).ToList();

                                    if (x3.Count <= max_features) features.AddRange(x3);
                                }
                            }

                            if (pse_aac_options.average_seq_position)
                            {
                                var x5 = average_seq_positions.Select(x => new feature_info()
                                {
                                    alphabet = alphabet_name,
                                    dimension = 1,
                                    category = $@"{category_prefix}_{nameof(average_seq_positions)}",
                                    source = source.ToString(),
                                    @group = $@"{group_prefix}_{sq.name}_{nameof(average_seq_positions)}_{alphabet_name}_{dist_name}",
                                    member = $@"{sq_index}_{x.name}",
                                    perspective = $@"default",
                                    feature_value = x.value
                                }).ToList();
                                if (x5.Count <= max_features) features.AddRange(x5);
                            }

                            if (pse_aac_options.average_dipeptide_distance)
                            {
                                var x6 = average_dipeptide_distance.Select(x => new feature_info()
                                {
                                    alphabet = alphabet_name,
                                    dimension = 1,
                                    category = $@"{category_prefix}_{nameof(average_dipeptide_distance)}",
                                    source = source.ToString(),
                                    @group = $@"{group_prefix}_{sq.name}_{nameof(average_dipeptide_distance)}_{alphabet_name}_{dist_name}",
                                    member = $@"{sq_index}_{x.name}",
                                    perspective = $@"default",
                                    feature_value = x.value
                                }).ToList();
                                if (x6.Count <= max_features) features.AddRange(x6);
                            }
                        }
                    }
                }
            }

            if (calculate_aa_or_ss_sequence_classification_data_template == null)
            {
                var template = features.Select(a => new feature_info(a) { source = "", feature_value = 0 }).ToList();
                calculate_aa_or_ss_sequence_classification_data_template = template;
            }

            return features;
        }


        public class feature_info_container
        {
            public List<feature_info> feautre_info_list;

 
            public static string serialise_json(feature_info_container feature_info_container)
            {    
                var json_settings = new JsonSerializerSettings() { PreserveReferencesHandling = PreserveReferencesHandling.All, };
                var serialised_serialise_json = JsonConvert.SerializeObject(feature_info_container, json_settings);
                return serialised_serialise_json;
            }

            public static feature_info_container deserialise(string serialized_json)
            {
                var json_settings = new JsonSerializerSettings() { PreserveReferencesHandling = PreserveReferencesHandling.All };
                feature_info_container feature_info_container = null;

                try
                {
                    feature_info_container = JsonConvert.DeserializeObject<feature_info_container>(serialized_json, json_settings);
                }
                catch (Exception e)
                {
                    feature_info_container = null;

                    Program.WriteLine(e.ToString());
                }

                return feature_info_container;
            }
        }

        public class feature_info
        {
            public string alphabet;
            public int dimension;
            public string category;
            public string source;
            public string group;
            public string member;
            public string perspective;
            public double feature_value;

            public feature_info()
            {

            }

            public feature_info(feature_info feature_info)
            {
                this.alphabet = feature_info.alphabet;
                this.dimension = feature_info.dimension;
                this.category = feature_info.category;
                this.source = feature_info.source;
                this.@group = feature_info.@group;
                this.member = feature_info.member;
                this.perspective = feature_info.perspective;
                this.feature_value = feature_info.feature_value;
            }

            public bool validate()
            {
                var ok = (!string.IsNullOrWhiteSpace(alphabet) &&
                    (dimension >= 0 && dimension <= 3) &&
                    !string.IsNullOrWhiteSpace(source) &&
                    !string.IsNullOrWhiteSpace(group) &&
                    !string.IsNullOrWhiteSpace(member) &&
                    !string.IsNullOrWhiteSpace(perspective) &&
                    (!double.IsNaN(feature_value) && !double.IsInfinity(feature_value)));

                return ok;
            }

            public override string ToString()
            {
                var data = new List<(string key, string value)>()
                {
                    (nameof(alphabet), alphabet),
                    (nameof(dimension), dimension.ToString()),
                    (nameof(category), category),
                    (nameof(source), source),
                    (nameof(group), group),
                    (nameof(member), member),
                    (nameof(perspective), perspective),
                    (nameof(feature_value), feature_value.ToString())
                };

                return string.Join(", ", data);
            }
        }

        //public static List<(instance_meta_data instance_meta_data, List<feature_info>)> encode_sequence_list(int class_id, protein_data_sources source, List<instance_meta_data> aa_sequences, int max_features)
        //{
        //    var pse_aac_options = new feature_calcs.pse_aac_options()
        //    {
        //        saac = true,
        //        oaac = true,
        //        motifs = true,
        //        dipeptides = true,
        //        average_seq_position = true,
        //        average_dipeptide_distance = true,
        //        dipeptides_binary = true,
        //        oaac_binary = true,
        //        saac_binary = true,
        //        motifs_binary = true,
        //    };
        //
        //    var features = aa_sequences.Select(aa_seq => (aa_seq, x: calculate_aa_or_ss_sequence_classification_data(source, "aa", "aa", aa_seq.subsequence_aa_seq, feature_calcs.seq_type.amino_acid_sequence, pse_aac_options, max_features))).ToList();
        //
        //    var class_id_feature = new feature_info()
        //    {
        //        alphabet = null,
        //        dimension = 0,
        //        category = nameof(class_id),
        //        source = nameof(class_id),
        //        group = nameof(class_id),
        //        member = nameof(class_id),
        //        perspective = nameof(class_id),
        //        feature_value = (double)class_id,
        //    };
        //
        //    for (var i = 0; i < features.Count; i++)
        //    {
        //        features[i].x.Insert(0, class_id_feature);
        //    }
        //
        //    //features = features.Select(a => (a.dimer_type, a.parallelism, a.symmetry_mode, a.pdb_id, a.chain_id, a.class_id, a.subsequence, a.x.OrderBy(b => b.source).ThenBy(b => b.group).ThenBy(b => b.member).ThenBy(b => b.perspective).ToList())).ToList();
        //
        //    return features;
        //}

        public static void copy_nosd(List<feature_info> features)
        {
            var copy_without_std_dev = true;

            if (copy_without_std_dev)
            {
                // find which groups contain dev_standard
                var groups_with_std_dev = features.GroupBy(a => (a.source, a.@group)).Where(a => a.Any(b => b.perspective == nameof(descriptive_stats.dev_standard))).SelectMany(a => a.ToList()).ToList();

                var copy_without_sd_perspectives = groups_with_std_dev.Where(a => a.perspective != nameof(descriptive_stats.dev_standard)).Select(a => new feature_info(a)
                {
                    @group = $"{a.@group}_nosd"
                }).ToList();

                features.AddRange(copy_without_sd_perspectives);
            }
        }

        public enum protein_data_sources
        {
            subsequence_1d,
            neighbourhood_1d,
            protein_1d,
            subsequence_3d,
            neighbourhood_3d,
            protein_3d
        }

        public static (instance_meta_data instance_meta_data, List<feature_info> feature_info) encode_subsequence_classification_data_row(
            subsequence_classification_data scd,
            int max_features,
            feature_types_1d feature_types_subsequence_1d,
            feature_types_1d feature_types_neighbourhood_1d,
            feature_types_1d feature_types_protein_1d,
            feature_types_3d feature_types_subsequence_3d,
            feature_types_3d feature_types_neighbourhood_3d,
            feature_types_3d feature_types_protein_3d)
        {

#if DEBUG
            if (Program.verbose_debug) Program.WriteLine($"{nameof(encode_subsequence_classification_data_row)}(sdc, max_features, feature_types_subsequence_1d, feature_types_neighbourhood_1d, feature_types_protein_1d, feature_types_subsequence_3d, feature_types_neighbourhood_3d, feature_types_protein_3d)");
#endif
            var check_num_features_consistency = true;

            feature_info class_id = null;



            List<feature_info> subsequence_1d_classification_data = null;
            List<feature_info> neighbourhood_1d_classification_data = null;
            List<feature_info> protein_1d_classification_data = null;

            List<feature_info> subsequence_3d_classification_data = null;
            List<feature_info> neighbourhood_3d_classification_data = null;
            List<feature_info> protein_3d_classification_data = null;

            class_id = calculate_class_id_classification_data(scd);


            var tasks = new List<Task>();

            if (feature_types_subsequence_1d != null)
            {
                var task = Task.Run(() =>
                {
                    subsequence_1d_classification_data = calculate_classification_data_1d(scd, scd.subsequence_master_atoms, protein_data_sources.subsequence_1d, max_features, feature_types_subsequence_1d);

                    subsequence_1d_classification_data = subsequence_1d_classification_data.GroupBy(a => (a.source, a.alphabet, a.category, a.dimension, a.@group)).Where(a => a.Count() <= max_features).SelectMany(a => a).ToList();

                    if (check_num_features_consistency)
                    {
                        var total_subsequence_1d_classification_data = subsequence_1d_classification_data?.Count ?? -1;
                        if (total_subsequence_1d_classification_data > -1 && subsequence_classification_data.total_subsequence_1d_classification_data > -1 && subsequence_classification_data.total_subsequence_1d_classification_data != total_subsequence_1d_classification_data) throw new Exception();
                        if (total_subsequence_1d_classification_data > -1) subsequence_classification_data.total_subsequence_1d_classification_data = total_subsequence_1d_classification_data;
                    }
                });
                tasks.Add(task);
            }

            if (feature_types_neighbourhood_1d != null)
            {
                var task = Task.Run(() =>
                {
                    scd.neighbourhood_1d = Atom.get_intramolecular_neighbourhood_1d(scd);

                    if (scd.neighbourhood_1d.subsequence_atoms.Count == 0)
                    {
                        Program.WriteLine("Warning: " + scd.pdb_id + scd.chain_id + " (class " + scd.class_id + ") has no 1d neighbourhood data");
                    }

                    neighbourhood_1d_classification_data = calculate_classification_data_1d(scd.neighbourhood_1d, scd.neighbourhood_1d.subsequence_master_atoms, protein_data_sources.neighbourhood_1d, max_features, feature_types_neighbourhood_1d);

                    neighbourhood_1d_classification_data = neighbourhood_1d_classification_data.GroupBy(a => (a.source, a.alphabet, a.category, a.dimension, a.@group)).Where(a => a.Count() <= max_features).SelectMany(a => a).ToList();

                    if (check_num_features_consistency)
                    {
                        var total_neighbourhood_1d_classification_data = neighbourhood_1d_classification_data?.Count ?? -1;
                        if (total_neighbourhood_1d_classification_data > -1 && subsequence_classification_data.total_neighbourhood_1d_classification_data > -1 && subsequence_classification_data.total_neighbourhood_1d_classification_data != total_neighbourhood_1d_classification_data) throw new Exception();
                        if (total_neighbourhood_1d_classification_data > -1) subsequence_classification_data.total_neighbourhood_1d_classification_data = total_neighbourhood_1d_classification_data;
                    }
                });
                tasks.Add(task);
            }

            if (feature_types_protein_1d != null)
            {
                var task = Task.Run(() =>
                {
                    scd.protein_1d = Atom.get_intramolecular_protein_1d(scd);

                    if (scd.protein_1d.subsequence_atoms.Count == 0)
                    {
                        Program.WriteLine("Warning: " + scd.pdb_id + scd.chain_id + " (class " + scd.class_id + ") has no 1d protein data");
                    }

                    protein_1d_classification_data = calculate_classification_data_1d(scd.protein_1d, scd.protein_1d.subsequence_master_atoms, protein_data_sources.protein_1d, max_features, feature_types_protein_1d);

                    protein_1d_classification_data = protein_1d_classification_data.GroupBy(a => (a.source, a.alphabet, a.category, a.dimension, a.@group)).Where(a => a.Count() <= max_features).SelectMany(a => a).ToList();

                    if (check_num_features_consistency)
                    {
                        var total_protein_1d_classification_data = protein_1d_classification_data?.Count ?? -1;
                        if (total_protein_1d_classification_data > -1 && subsequence_classification_data.total_protein_1d_classification_data > -1 && subsequence_classification_data.total_protein_1d_classification_data != total_protein_1d_classification_data) throw new Exception();
                        if (total_protein_1d_classification_data > -1) subsequence_classification_data.total_protein_1d_classification_data = total_protein_1d_classification_data;
                    }
                });
                tasks.Add(task);
            }

            if (feature_types_subsequence_3d != null)
            {
                var task = Task.Run(() =>
                {

                    subsequence_3d_classification_data = calculate_classification_data_3d(scd, scd.subsequence_master_atoms, protein_data_sources.subsequence_3d, max_features, feature_types_subsequence_3d);

                    subsequence_3d_classification_data = subsequence_3d_classification_data.GroupBy(a => (a.source, a.alphabet, a.category, a.dimension, a.@group)).Where(a => a.Count() <= max_features).SelectMany(a => a).ToList();

                    if (check_num_features_consistency)
                    {
                        var total_subsequence_3d_classification_data = subsequence_3d_classification_data?.Count ?? -1;
                        if (total_subsequence_3d_classification_data > -1 && subsequence_classification_data.total_subsequence_3d_classification_data > -1 && subsequence_classification_data.total_subsequence_3d_classification_data != total_subsequence_3d_classification_data) throw new Exception();
                        if (total_subsequence_3d_classification_data > -1) subsequence_classification_data.total_subsequence_3d_classification_data = total_subsequence_3d_classification_data;
                    }
                });
                tasks.Add(task);

            }

            if (feature_types_neighbourhood_3d != null)
            {
                var task = Task.Run(() =>
                {
                    scd.neighbourhood_3d = Atom.get_intramolecular_neighbourhood_3d(scd);

                    if (scd.neighbourhood_3d.subsequence_atoms.Count == 0)
                    {
                        Program.WriteLine("Warning: " + scd.pdb_id + scd.chain_id + " (class " + scd.class_id + ") has no 3d neighbourhood data");
                    }

                    neighbourhood_3d_classification_data = calculate_classification_data_3d(scd.neighbourhood_3d, scd.neighbourhood_3d.subsequence_master_atoms, protein_data_sources.neighbourhood_3d, max_features, feature_types_neighbourhood_3d);

                    neighbourhood_3d_classification_data = neighbourhood_3d_classification_data.GroupBy(a => (a.source, a.alphabet, a.category, a.dimension, a.@group)).Where(a => a.Count() <= max_features).SelectMany(a => a).ToList();

                    if (check_num_features_consistency)
                    {
                        var total_neighbourhood_3d_classification_data = neighbourhood_3d_classification_data?.Count ?? -1;
                        if (total_neighbourhood_3d_classification_data > -1 && subsequence_classification_data.total_neighbourhood_3d_classification_data > -1 && subsequence_classification_data.total_neighbourhood_3d_classification_data != total_neighbourhood_3d_classification_data) throw new Exception();
                        if (total_neighbourhood_3d_classification_data > -1) subsequence_classification_data.total_neighbourhood_3d_classification_data = total_neighbourhood_3d_classification_data;
                    }
                });
                tasks.Add(task);

            }

            if (feature_types_protein_3d != null)
            {
                var task = Task.Run(() =>
                {
                    scd.protein_3d = Atom.get_intramolecular_protein_3d(scd);

                    if (scd.protein_3d.subsequence_atoms.Count == 0)
                    {
                        Program.WriteLine("Warning: " + scd.pdb_id + scd.chain_id + " (class " + scd.class_id + ") has no 3d protein data");
                    }

                    protein_3d_classification_data = calculate_classification_data_3d(scd.protein_3d, scd.protein_3d.subsequence_master_atoms, protein_data_sources.protein_3d, max_features, feature_types_protein_3d);

                    protein_3d_classification_data = protein_3d_classification_data.GroupBy(a => (a.source, a.alphabet, a.category, a.dimension, a.@group)).Where(a => a.Count() <= max_features).SelectMany(a => a).ToList();

                    if (check_num_features_consistency)
                    {
                        var total_protein_3d_classification_data = protein_3d_classification_data?.Count ?? -1;
                        if (total_protein_3d_classification_data > -1 && subsequence_classification_data.total_protein_3d_classification_data > -1 && subsequence_classification_data.total_protein_3d_classification_data != total_protein_3d_classification_data) throw new Exception();
                        if (total_protein_3d_classification_data > -1) subsequence_classification_data.total_protein_3d_classification_data = total_protein_3d_classification_data;
                    }
                });
                tasks.Add(task);
            }

            Task.WaitAll(tasks.ToArray<Task>());

            List<feature_info> features = new List<feature_info>();
            features.Add(class_id);
            if (subsequence_1d_classification_data != null && subsequence_1d_classification_data.Count > 0) features.AddRange(subsequence_1d_classification_data);
            if (neighbourhood_1d_classification_data != null && neighbourhood_1d_classification_data.Count > 0) features.AddRange(neighbourhood_1d_classification_data);
            if (protein_1d_classification_data != null && protein_1d_classification_data.Count > 0) features.AddRange(protein_1d_classification_data);
            if (subsequence_3d_classification_data != null && subsequence_3d_classification_data.Count > 0) features.AddRange(subsequence_3d_classification_data);
            if (neighbourhood_3d_classification_data != null && neighbourhood_3d_classification_data.Count > 0) features.AddRange(neighbourhood_3d_classification_data);
            if (protein_3d_classification_data != null && protein_3d_classification_data.Count > 0) features.AddRange(protein_3d_classification_data);

            var nosd = false;
            if (nosd)
            {
                copy_nosd(features);
            }

            var protein_seq = string.Join("", scd.pdb_chain_master_atoms.Select(a => a.amino_acid).ToList());
            var instance_meta_data = new instance_meta_data()
            {
                protein_pdb_id = scd.pdb_id,
                protein_chain_id = scd.chain_id,
                protein_dimer_type = scd.dimer_type,
                protein_aa_seq = protein_seq,

                subsequence_class_id = scd.class_id,
                subsequence_parallelism = scd.parallelism,
                subsequence_symmetry_mode = scd.symmetry_mode,
                subsequence_aa_seq = scd.aa_subsequence,
                subsequence_res_ids = scd.res_ids,
                subsequence_dssp_monomer = scd.dssp_monomer_subsequence,
                subsequence_dssp_multimer = scd.dssp_multimer_subsequence,
                subsequence_dssp3_monomer = Atom.secondary_structure_state_reduction(scd.dssp_monomer_subsequence),
                subsequence_dssp3_multimer = Atom.secondary_structure_state_reduction(scd.dssp_multimer_subsequence),
                ss_predictions = Atom.get_dssp_and_mpsa_subsequences(scd.pdb_id, scd.chain_id, scd.subsequence_master_atoms),
            };

            return (instance_meta_data, features);
        }

    }

}
