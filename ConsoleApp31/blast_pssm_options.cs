using System;

namespace Coil_DHC_Dataset_Maker
{

    public class blast_pssm_options
    {
        public bool make_unsplit_sequence = true;
        public bool make_split_sequence = true;

        public bool make_standard_encoding = true;
        public bool make_distance_transform = false;

        public bool normalise_none = false;
        public bool normalise_whole_pssm = true;
        public bool normalise_subsequence = false;
        public bool normalise_encoded_parts = false;
        public bool normalise_encoded_vector = false;

        public bool encode_standard_vector = true;
        public bool encoded_interval_vector = false;

        // separate groups with the following statistics:
        public bool encode_min = false;
        public bool encode_max = false;
        public bool encode_mean = true;
        public bool encode_mean_sd = false;
        public bool encode_median = false;
        public bool encode_mode = false;
        public bool encode_range = false;

        // encoding size (when using standard AA alphabet) - e.g. 20x20 = 400 could be 4x4 = 16 for reduced alphabet
        public bool size_1 = true;
        public bool size_20 = true;
        public bool size_210 = true;
        public bool size_400 = true;

        public bool db_nr_local_1e_4 = false; // 0.0001
        public bool db_nr_local_def = true; 
        public bool db_nr_remote_1e_4 = false;

        public bool db_sp_local_1e_4 = false;
        public bool db_sp_local_def = true;
        public bool db_sp_remote_1e_4 = false;

        public bool db_ur90_local_def = true;


        public void set_all(bool value)
        {
            make_unsplit_sequence = value;
            make_split_sequence = value;

            make_standard_encoding = value;
            make_distance_transform = value;

            normalise_none = value;
            normalise_whole_pssm = value;
            normalise_subsequence = value;
            normalise_encoded_parts = value;
            normalise_encoded_vector = value;

            encode_standard_vector = value;
            encoded_interval_vector = value;


            encode_min = value;
            encode_max = value;
            encode_mean = value;
            encode_mean_sd = value;
            encode_median = value;
            encode_mode = value;
            encode_range = value;


            size_1 = value;
            size_20 = value;
            size_210 = value;
            size_400 = value;

            db_nr_local_1e_4 = value;
            db_nr_local_def = value;
            db_nr_remote_1e_4 = value;

            db_sp_local_1e_4 = value;
            db_sp_local_def = value;
            db_sp_remote_1e_4 = value;

            db_ur90_local_def = value;
        }
    }

}