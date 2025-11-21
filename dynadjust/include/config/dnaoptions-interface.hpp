//============================================================================
// Name         : dnaoptions-interface.hpp
// Author       : Roger Fraser
// Contributors :
// Version      : 1.00
// Copyright    : Copyright 2017 Geoscience Australia
//
//                Licensed under the Apache License, Version 2.0 (the "License");
//                you may not use this file except in compliance with the License.
//                You may obtain a copy of the License at
//               
//                http ://www.apache.org/licenses/LICENSE-2.0
//               
//                Unless required by applicable law or agreed to in writing, software
//                distributed under the License is distributed on an "AS IS" BASIS,
//                WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//                See the License for the specific language governing permissions and
//                limitations under the License.
//
// Description  : DynAdjust interface options include file
//============================================================================

#ifndef DNAOPTIONS_INTERFACE_HPP
#define DNAOPTIONS_INTERFACE_HPP

#if defined(_MSC_VER)
	#if defined(LIST_INCLUDES_ON_BUILD) 
		#pragma message("  " __FILE__) 
	#endif
#endif

const char* const ALL_MODULE_STDOPT = "Standard options";
const char* const ALL_MODULE_OUTPUT = "Output options";
const char* const ALL_MODULE_EXPORT = "Export options";
const char* const ALL_MODULE_GENERIC = "Generic options";

const char* const IMPORT_MODULE_FRAME = "Reference frame options";
const char* const IMPORT_MODULE_SCREEN = "Data screening options";
const char* const IMPORT_MODULE_GNSS_VAR = "GNSS variance matrix scaling options";
const char* const IMPORT_MODULE_SIMULATE = "Network simulation options";

const char* const ADJUST_MODULE_MODE = "Adjustment mode options";
const char* const ADJUST_MODULE_PHASED = "Phased adjustment options";
const char* const ADJUST_MODULE_CONFIG = "Configuration options";
const char* const ADJUST_MODULE_STAGE = "Stage adjustment options";

const char* const SEGMENT_MODULE_CONFIG = "Configuration options";

const char* const GEOID_MODULE_NTV2 = "NTv2 creation options";
const char* const GEOID_MODULE_INTERPOLATE = "Interpolation options";
const char* const GEOID_MODULE_INTERACTIVE = "Interactive interpolation options";
const char* const GEOID_MODULE_FILE = "Text file interpolation options";

const char* const REFTRAN_MODULE_TRANSFORM = "Transformation options";

const char* const PLOT_MODULE_CONFIG = "Data configuration options";
const char* const PLOT_MODULE_TITLEBLOCK = "Title block options";
const char* const PLOT_MODULE_MAP = "Mapping options";

const char* const IMPORT_FILE = "stn-msr-file";
const char* const IMPORT_FILE_F = "stn-msr-file,f";
const char* const NETWORK_NAME = "network-name";
const char* const NETWORK_NAME_N = "network-name,n";
const char* const INPUT_FOLDER = "input-folder";
const char* const INPUT_FOLDER_I = "input-folder,i";
const char* const OUTPUT_FOLDER = "output-folder";
const char* const OUTPUT_FOLDER_O = "output-folder,o";
const char* const DYNADJUST_LOG_FILE = "dynadjust-log-file";
const char* const PROJECT_FILE = "project-file";
const char* const PROJECT_FILE_P = "project-file,p";
const char* const BIN_STN_FILE = "binary-stn-file";
const char* const BIN_STN_FILE_S = "binary-stn-file,s";
const char* const BIN_MSR_FILE = "binary-msr-file";
const char* const BIN_MSR_FILE_M = "binary-msr-file,m";
const char* const CONVERT_BST_HT = "convert-stn-hts";
const char* const HELP_MODULE = "help-module";
const char* const HELP = "help";
const char* const HELP_H = "help,h";
const char* const VERSION = "version";
const char* const VERSION_V = "version,v";
const char* const QUIET = "quiet";
const char* const VERBOSE = "verbose-level";
const char* const COUT = "cout";
const char* const CHECK_INPUT_COORD = "check-input-coord";

const char* const RUN_IMPORT = "import";
const char* const RUN_GEOID = "geoid";
const char* const RUN_REFTRAN = "reftran";
const char* const RUN_SEGMENT = "segment";
const char* const RUN_ADJUST = "adjust";

const char* const VARIABLES = "variables";

const char* const ADJUSTMENT_MODE = "adjustment-mode";
const char* const MODE_PHASED = "phased-adjustment";
const char* const STAGED_ADJUSTMENT = "staged-adjustment";
const char* const MODE_PHASED_BLOCK1 = "block1-phased";
const char* const MODE_SIMULTANEOUS = "simultaneous-adjustment";
const char* const MODE_SIMULATION = "simulation";
const char* const MODE_ADJ_REPORT = "report-results";
const char* const MODE_PHASED_MT = "multi-thread";

const char* const COMMENTS = "comments";
const char* const CONF_INTERVAL = "conf-interval";
const char* const ITERATION_THRESHOLD = "iteration-threshold";
const char* const MAX_ITERATIONS = "max-iterations";
const char* const FREE_STN_SD = "free-stn-sd";
const char* const FIXED_STN_SD = "fixed-stn-sd";
const char* const STN_CONSTRAINTS = "constraints";
const char* const OUTPUT_ADJ_STN_ITER =  "output-iter-adj-stn";
const char* const OUTPUT_ADJ_STAT_ITER = "output-iter-adj-stat";
const char* const OUTPUT_ADJ_MSR_ITER =  "output-iter-adj-msr";
const char* const OUTPUT_CMP_MSR_ITER =  "output-iter-cmp-msr";
const char* const OUTPUT_ADJ_MSR = "output-adj-msr";
const char* const OUTPUT_ADJ_GNSS_UNITS = "output-adj-gnss-units";
const char* const OUTPUT_ADJ_MSR_TSTAT = "output-tstat-adj-msr";
const char* const OUTPUT_STN_COR_FILE = "output-corrections-file";
const char* const OUTPUT_STN_CORR = "stn-corrections";
const char* const HZ_CORR_THRESHOLD = "hz-corr-threshold";
const char* const VT_CORR_THRESHOLD = "vt-corr-threshold";
const char* const OUTPUT_POS_UNCERTAINTY = "output-pos-uncertainty";
const char* const OUTPUT_APU_CORRELATIONS = "output-all-covariances";
const char* const OUTPUT_APU_UNITS = "output-apu-vcv-units";
const char* const OUTPUT_MSR_TO_STN = "output-msr-to-stn";

const char* const TYPE_B_GLOBAL = "type-b-sd-global";
const char* const TYPE_B_FILE = "type-b-sd-file";

const char* const OUTPUT_ADJ_MSR_SORTBY = "sort-adj-msr-field";
const char* const OUTPUT_MSR_TO_STN_SORTBY = "sort-msr-to-stn-field";

const char* const OUTPUT_ADJ_MSR_DBID = "output-database-ids";

const char* const OUTPUT_IGNORED_MSRS = "output-ignored-msrs";

const char* const OUTPUT_ADJ_STN_BLOCKS = "output-stn-blocks";
const char* const OUTPUT_ADJ_MSR_BLOCKS = "output-msr-blocks";
const char* const OUTPUT_ADJ_STN_SORT_ORDER = "sort-stn-orig-order";
const char* const OUTPUT_STN_COORD_TYPES = "stn-coord-types";
const char* const OUTPUT_PRECISION_SECONDS_MSR = "precision-msr-angular";
const char* const OUTPUT_PRECISION_METRES_MSR = "precision-msr-linear";
const char* const OUTPUT_PRECISION_SECONDS_STN = "precision-stn-angular";
const char* const OUTPUT_PRECISION_METRES_STN = "precision-stn-linear";
const char* const OUTPUT_ANGULAR_TYPE_MSR = "angular-msr-type";
const char* const OUTPUT_DMS_FORMAT_MSR = "dms-msr-format";
const char* const OUTPUT_ANGULAR_TYPE_STN = "angular-stn-type";

//const char* const MVAR_INVERSE_METHOD = "msr-inverse-method";
const char* const LSQ_INVERSE_METHOD = "inversion-method";
const char* const SCALE_NORMAL_UNITY = "scale-normals-to-unity";
const char* const PURGE_STAGE_FILES = "purge-stage-files";
const char* const RECREATE_STAGE_FILES = "create-stage-files";
const char* const UPDATE_ORIGINAL_STN_FILE = "update-orig-stn-file";

const char* const SEG_MIN_INNER_STNS = "min-inner-stns";
const char* const SEG_THRESHOLD_STNS = "max-block-stns";
const char* const SEG_STARTING_STN = "starting-stns";
const char* const SEG_SORT_STNS_MSR = "sort-block-stns";
const char* const SEG_VIEW_ON_SEGMENT = "view-block-on-segment";
const char* const SEG_PRINT_DEBUG = "print-segment-debug";
const char* const SEG_SHOW_SUMMARY = "show-segment-summary";
const char* const SEG_DISP_BLK_NET = "display-block-network";
const char* const SEG_FORCE_CONTIGUOUS = "contiguous-blocks";
const char* const SEG_SEARCH_LEVEL = "search-level";

const char* const GEOID_PATH = "geoid-file";
const char* const INTERPOLATE_ALWAYS = "interpolate-heights-always";

const char* const REFERENCE_FRAME = "reference-frame";
const char* const REFERENCE_FRAME_R = "reference-frame,r";
const char* const EPOCH = "epoch";
const char* const EPOCH_E = "epoch,e";
const char* const OVERRIDE_INPUT_FRAME = "override-input-ref-frame";
const char* const TECTONIC_PLATE_BDY_FILE = "plate-boundary-file";
const char* const TECTONIC_PLATE_BDY_FILE_B = "plate-boundary-file,b";
const char* const TECTONIC_PLATE_POLE_FILE = "plate-pole-file";
const char* const TECTONIC_PLATE_POLE_FILE_M = "plate-pole-file,m";
const char* const TECTONIC_PLATE_MODEL_OPTION = "plate-model-option";

const char* const VSCALE = "v-scale";
const char* const PSCALE = "p-scale";
const char* const LSCALE = "l-scale";
const char* const HSCALE = "h-scale";
const char* const SCALAR_FILE = "baseline-scalar-file";

const char* const EXPORT_XML_FILES = "export-xml-files";
const char* const EXPORT_XML_STN_FILE = "export-xml-stn-file";
const char* const EXPORT_XML_MSR_FILE = "export-xml-msr-file";
const char* const EXPORT_FROM_BINARY = "export-from-binary";
const char* const EXPORT_SINGLE_XML_FILE = "single-xml-file";
const char* const PREFER_X_MSR_AS_G = "prefer-single-x-as-g";
const char* const EXPORT_DNA_FILES = "export-dna-files";
const char* const EXPORT_DNA_STN_FILE = "export-dna-stn-file";
const char* const EXPORT_DNA_MSR_FILE = "export-dna-msr-file";
const char* const EXPORT_GEO_FILE = "export-dna-geo-file";
const char* const EXPORT_NTV2_ASCII_FILE = "export-ntv2-asc-file";
const char* const EXPORT_NTV2_BINARY_FILE = "export-ntv2-gsb-file";
const char* const EXPORT_SNX_FILE = "export-sinex-file";
const char* const IMPORT_GEO_FILE = "geo-file";
const char* const IMPORT_GEO_FILE_G = "geo-file,g";
const char* const EXPORT_ASL_FILE = "export-asl-file";
const char* const EXPORT_AML_FILE = "export-aml-file";
const char* const EXPORT_MAP_FILE = "export-map-file";
const char* const EXPORT_DISCONT_FILE = "export-discont-file";
const char* const SIMULATE_MSR_FILE = "simulate-msr-file";
const char* const STATION_RENAMING_FILE = "stn-renaming-file";
const char* const STATION_DISCONTINUITY_FILE = "discontinuity-file";
const char* const TEST_INTEGRITY = "test-integrity";
const char* const VERIFY_COORDS = "verify-coordinates";
const char* const INCLUDE_MSRS = "include-msr-types";
const char* const EXCLUDE_MSRS = "exclude-msr-types";
const char* const BOUNDING_BOX = "bounding-box";
const char* const GET_MSRS_TRANSCENDING_BOX = "get-msrs-transcending-box";
const char* const INCLUDE_STN_ASSOC_MSRS =    "include-stns-assoc-msrs";
const char* const EXCLUDE_STN_ASSOC_MSRS = "exclude-stns-assoc-msrs";
const char* const SPLIT_CLUSTERS = "split-gnss-cluster-msrs";
const char* const TEST_NEARBY_STNS = "search-nearby-stn";
const char* const TEST_NEARBY_STN_DIST = "nearby-stn-buffer";
const char* const TEST_SIMILAR_MSRS = "search-similar-msr";
const char* const TEST_SIMILAR_GNSS_MSRS = "search-similar-gnss-msr";
const char* const IGNORE_SIMILAR_MSRS = "ignore-similar-msr";
const char* const REMOVE_IGNORED_MSRS = "remove-ignored-msr";
const char* const FLAG_UNUSED_STNS = "flag-unused-stations";
const char* const IMPORT_SEG_BLOCK = "import-block-stn-msr";

const char* const PLOT_MSRS = "plot-msr-types";
const char* const PLOT_MSRS_IGNORED = "plot-ignored-msrs";
const char* const PLOT_STN_LABELS = "label-stations";
const char* const PLOT_ALT_NAME = "alternate-name";			
const char* const PLOT_CONSTRAINT_LABELS = "label-constraints";
const char* const PLOT_CORRECTION_ARROWS = "correction-arrows";
const char* const PLOT_ERROR_ELLIPSES = "error-ellipses";
const char* const PLOT_POSITIONAL_UNCERTAINTY = "positional-uncertainty";
const char* const CORRECTION_SCALE = "scale-arrows";
const char* const PU_ELLIPSE_SCALE = "scale-ellipse-circles";
const char* const PLOT_CORRECTION_LABELS = "label-corrections";
const char* const COMPUTE_CORRECTIONS = "compute-corrections";
const char* const OMIT_MEASUREMENTS = "omit-measurements";
const char* const PROJECTION = "map-projection";
const char* const LABEL_FONT_SIZE = "label-font-size";
const char* const MSR_LINE_WIDTH = "msr-line-width";
const char* const PLOT_BLOCKS = "phased-block-view";
const char* const PLOT_CENTRE_LAT = "centre-latitude";
const char* const PLOT_CENTRE_LON = "centre-longitude";
const char* const PLOT_CENTRE_STATION = "centre-station";
const char* const PLOT_AREA_RADIUS = "area-radius";
const char* const OMIT_TITLE_BLOCK = "omit-title-block";
const char* const DONT_CREATE_PDF = "supress-pdf-creation";
const char* const GRAPH_SEGMENTATION_STNS = "graph-stn-blocks";
const char* const GRAPH_SEGMENTATION_MSRS = "graph-msr-blocks";
const char* const BLOCK_NUMBER = "block-number";
const char* const TITLEBLOCK_TITLE = "title";
const char* const TITLEBLOCK_NAME = "org-unit-name";
const char* const TITLEBLOCK_SUBNAME = "org-subunit-name";
const char* const KEEP_FILES = "keep-gen-files";
const char* const EXPORT_PNG = "export-png";

const char* const INTERACTIVE = "interactive";
const char* const INTERACTIVE_E = "interactive,e";
// const char* const FILE_INTERPOLATE = "file-interpolate";
// const char* const FILE_INTERPOLATE_F = "file-interpolate,f";
const char* const CREATE_NTV2 = "create-ntv2";
const char* const CREATE_NTV2_C = "create-ntv2,c";
const char* const SUMMARY = "summary";
const char* const SUMMARY_U = "summary,u";
const char* const NTV2_FILEPATH = "ntv2-file";
const char* const NTV2_FILEPATH_G = "ntv2-file,g";
const char* const DAT_FILEPATH = "dat-file";
const char* const DAT_FILEPATH_D = "dat-file,d";
const char* const LATITUDE = "latitude";
const char* const LATITUDE_P = "latitude,p";
const char* const LONGITUDE = "longitude";
const char* const LONGITUDE_L = "longitude,l";
const char* const INPUT_FILE = "text-file";
const char* const INPUT_FILE_T = "text-file,t";
const char* const DIRECTION = "direction";
const char* const DIRECTION_R = "direction,r";
const char* const METHOD = "interpolation-method";
const char* const METHOD_M = "interpolation-method,m";
const char* const DDEG_FORMAT = "decimal-degrees";

const char* const NTV2_GS_TYPE = "grid-shift-type";
const char* const NTV2_VERSION = "grid-version";
const char* const NTV2_SYSTEM_F = "system-from";
const char* const NTV2_SYSTEM_T = "system-to";
const char* const NTV2_MAJOR_F = "semi-major-from";
const char* const NTV2_MAJOR_T = "semi-major-to";
const char* const NTV2_MINOR_F = "semi-minor-from";
const char* const NTV2_MINOR_T = "semi-minor-to";
const char* const NTV2_SUB_NAME = "sub-grid-name";
const char* const NTV2_CREATED = "creation-date";
const char* const NTV2_UPDATED = "update-date";

const char* const MAP_FILE_OPT = "map-file";
const char* const ASL_FILE = "asl-file";
const char* const AML_FILE = "aml-file";
const char* const DST_FILE = "dst-file";
const char* const DMS_FILE = "dms-file";
const char* const XMLOUT_FILE = "xml-file";
const char* const SEG_FILE = "seg-file";
const char* const NET_FILE = "net-file";
const char* const XYZ_FILE = "xyz-file";
const char* const ADJ_FILE = "adj-file";
const char* const SNX_FILE = "snx-file";

const char* const section_variables = "#variables";
const char* const section_general = "#general";
const char* const section_import = "#import";
const char* const section_reftran = "#reftran";
const char* const section_geoid = "#geoid";
const char* const section_segment = "#segment";
const char* const section_adjust = "#adjust";
const char* const section_output = "#output";
const char* const section_plot = "#plot";
const char* const section_display = "#display";

#endif  // DNAOPTIONS_INTERFACE_HPP
