interface
   function initial_file_get_id_codon(tag) result(tag_result) bind(c, name="initial_file_get_id_codon")
     use iso_c_binding, only: c_int64_t
     integer(c_int64_t), value :: tag
     integer(c_int64_t) :: tag_result
   end function initial_file_get_id_codon
   function topo_file_get_id_codon(tag) result(tag_result) bind(c, name="topo_file_get_id_codon")
     use iso_c_binding, only: c_int64_t
     integer(c_int64_t), value :: tag
     integer(c_int64_t) :: tag_result
   end function topo_file_get_id_codon
   function getname_codon(tag) result(tag_out) bind(c, name="getname_codon")
     use iso_c_binding, only: c_int64_t
     integer(c_int64_t), value :: tag
     integer(c_int64_t) :: tag_out
   end function getname_codon
   function sec2hms_codon(tag) result(tag_out) bind(c, name="sec2hms_codon")
     use iso_c_binding, only: c_int64_t
     integer(c_int64_t), value :: tag
     integer(c_int64_t) :: tag_out
   end function sec2hms_codon
   function write_hist_coord_att_codon(tag) result(tag_out) bind(c, name="write_hist_coord_att_codon")
     use iso_c_binding, only: c_int64_t
     integer(c_int64_t), value :: tag
     integer(c_int64_t) :: tag_out
   end function write_hist_coord_att_codon
   function read_initial_codon(tag) result(tag_out) bind(c, name="read_initial_codon")
     use iso_c_binding, only: c_int64_t
     integer(c_int64_t), value :: tag
     integer(c_int64_t) :: tag_out
   end function read_initial_codon
   function preset_codon(tag) result(tag_out) bind(c, name="preset_codon")
     use iso_c_binding, only: c_int64_t
     integer(c_int64_t), value :: tag
     integer(c_int64_t) :: tag_out
   end function preset_codon
   function cam_initfiles_close_codon(tag) result(tag_out) bind(c, name="cam_initfiles_close_codon")
     use iso_c_binding, only: c_int64_t
     integer(c_int64_t), value :: tag
     integer(c_int64_t) :: tag_out
   end function cam_initfiles_close_codon
   function cam_initfiles_open_codon(tag) result(tag_out) bind(c, name="cam_initfiles_open_codon")
     use iso_c_binding, only: c_int64_t
     integer(c_int64_t), value :: tag
     integer(c_int64_t) :: tag_out
   end function cam_initfiles_open_codon
   function restart_setopts_codon(tag) result(tag_out) bind(c, name="restart_setopts_codon")
     use iso_c_binding, only: c_int64_t
     integer(c_int64_t), value :: tag
     integer(c_int64_t) :: tag_out
   end function restart_setopts_codon
   function write_inithist_codon(tag) result(tag_out) bind(c, name="write_inithist_codon")
     use iso_c_binding, only: c_int64_t
     integer(c_int64_t), value :: tag
     integer(c_int64_t) :: tag_out
   end function write_inithist_codon
   function lookup_hist_coord_indices_codon(tag) result(tag_out) bind(c, name="lookup_hist_coord_indices_codon")
     use iso_c_binding, only: c_int64_t
     integer(c_int64_t), value :: tag
     integer(c_int64_t) :: tag_out
   end function lookup_hist_coord_indices_codon
   function check_var_codon(tag) result(tag_out) bind(c, name="check_var_codon")
     use iso_c_binding, only: c_int64_t
     integer(c_int64_t), value :: tag
     integer(c_int64_t) :: tag_out
   end function check_var_codon
   function scam_default_opts_codon(tag) result(tag_out) bind(c, name="scam_default_opts_codon")
     use iso_c_binding, only: c_int64_t
     integer(c_int64_t), value :: tag
     integer(c_int64_t) :: tag_out
   end function scam_default_opts_codon
   function formula_terms_copy_codon(tag) result(tag_out) bind(c, name="formula_terms_copy_codon")
     use iso_c_binding, only: c_int64_t
     integer(c_int64_t), value :: tag
     integer(c_int64_t) :: tag_out
   end function formula_terms_copy_codon
   function get_masterlist_indx_codon(tag) result(tag_out) bind(c, name="get_masterlist_indx_codon")
     use iso_c_binding, only: c_int64_t
     integer(c_int64_t), value :: tag
     integer(c_int64_t) :: tag_out
   end function get_masterlist_indx_codon
   function bld_htapefld_indices_codon(tag) result(tag_out) bind(c, name="bld_htapefld_indices_codon")
     use iso_c_binding, only: c_int64_t
     integer(c_int64_t), value :: tag
     integer(c_int64_t) :: tag_out
   end function bld_htapefld_indices_codon
   function inifld_codon(tag) result(tag_out) bind(c, name="inifld_codon")
     use iso_c_binding, only: c_int64_t
     integer(c_int64_t), value :: tag
     integer(c_int64_t) :: tag_out
   end function inifld_codon
   function setup_interpolation_and_define_vector_compliments_codon(tag) result(tag_out) &
        bind(c, name="setup_interpolation_and_define_vector_compliments_codon")
     use iso_c_binding, only: c_int64_t
     integer(c_int64_t), value :: tag
     integer(c_int64_t) :: tag_out
   end function setup_interpolation_and_define_vector_compliments_codon
   function add_hist_coord_regonly_codon(tag) result(tag_out) bind(c, name="add_hist_coord_regonly_codon")
     use iso_c_binding, only: c_int64_t
     integer(c_int64_t), value :: tag
     integer(c_int64_t) :: tag_out
   end function add_hist_coord_regonly_codon
   function add_default_codon(tag) result(tag_out) bind(c, name="add_default_codon")
     use iso_c_binding, only: c_int64_t
     integer(c_int64_t), value :: tag
     integer(c_int64_t) :: tag_out
   end function add_default_codon
   function init_restart_history_codon(tag) result(tag_out) bind(c, name="init_restart_history_codon")
     use iso_c_binding, only: c_int64_t
     integer(c_int64_t), value :: tag
     integer(c_int64_t) :: tag_out
   end function init_restart_history_codon
end interface
