#!/usr/bin/env bash

# patch setup script with github version (conda version missing --env_name argument)

patch $CONDA_PREFIX/opt/aquamis/scripts/aquamis_setup.sh <<'EOF'
--- aquamis_setup-old.sh	2023-09-19 21:57:10.000000000 +0200
+++ aquamis_setup-new.sh	2023-09-26 09:30:14.000000000 +0200
@@ -30,6 +30,7 @@
 For more information, please visit https://gitlab.com/bfr_bioinformatics/AQUAMIS
 
 ${bold}Options:${normal}
+  -e, --env_name             Override conda environment name derived from envs/aquamis.yaml (default: aquamis)
   -m, --mamba                Install the latest version of 'mamba' to the Conda base environment and
                              create the AQUAMIS environment from the git repository recipe
   -b, --busco                Download Augustus and BUSCO databases to <AQUAMIS>/download and extract them in the conda environment
@@ -38,6 +39,7 @@
   -s, --status               Show installation status
   -f, --force                Force overwrite for downloads in <AQUAMIS>/download
   -k, --keep_downloads       Do not remove downloads after extraction
+  -a, --auto                 Do not ask for interactive confirmation
   -v, --verbose              Print script debug info
   -h, --help                 Show this help
 
@@ -54,7 +56,7 @@
   local args=("$@")
 
   ## default values of variables
-  conda_recipe="aquamis.yaml"
+  conda_recipe="$repo_path/envs/aquamis.yaml"
 
   ## default values of switches
   conda_status=false
@@ -68,6 +70,7 @@
   arg_status=false
   arg_force=false
   arg_keep_dl=false
+  interactive=true
 
   while :; do
     case "${1-}" in
@@ -76,6 +79,10 @@
         set -x
         shift
         ;;
+      -e | --env_name)
+        conda_recipe_env_name="$2"
+        shift 2
+        ;;
       -m | --mamba)
         arg_mamba=true
         shift
@@ -104,6 +111,10 @@
         arg_keep_dl=true
         shift
         ;;
+      -a | --auto)
+        interactive=false
+        shift
+        ;;
       -?*) die "Unknown option: $1. For the manual, type ${blue}bash $script_real --help${normal}" ;;
       *) break ;;
     esac
@@ -114,6 +125,7 @@
 
   ## checks
   [[ $(basename $script_path) != "scripts" ]] && die "This setup script does not reside in its original installation directory. This is a requirement for proper execution. Aborting..."
+  [[ -f $conda_recipe ]] || die "The conda recipe does not exist: $conda_recipe"
 
   return 0
 }
@@ -128,7 +140,10 @@
 show_info() {
   logheader "Parameters:"
   echo "
-repo_path:|$repo_path
+Repository:|$repo_path
+Conda Recipe:|$conda_recipe
+Conda Recipe Env Name:|$conda_recipe_env_name
+Keep Downloads:|$arg_keep_dl
 " | column -t -s "|"
 }
 
@@ -137,8 +152,8 @@
 download_file() {
   # Init
   download_success=false
-  download_status=false
-  download_hash=''
+  local download_status=false
+  local download_hash=''
   
   # Download
   if [[ -f $local_archive ]] && [[ $arg_force == false ]]; then
@@ -167,10 +182,10 @@
   [[ "$mamba_status" == true ]] && logwarn "Skipping Mamba installation. Mamba is already detected."
 
   # Create Conda Environment for AQUAMIS
-  logentry "Creating the AQUAMIS environment with the Conda recipe: $repo_path/envs/${conda_recipe}"
+  logentry "Creating the AQUAMIS environment with the Conda recipe: $conda_recipe"
   set +eu  # workaround: see https://github.com/conda/conda/issues/8186
-  [[ -z "${conda_env_path}" ]] && { mamba env create -f $repo_path/envs/${conda_recipe} || true; }
-  [[ -n "${conda_env_path}" && "$arg_force" == true && ! -d "${conda_base}/envs/${conda_env_name}" ]] && { mamba env create -f $repo_path/envs/${conda_recipe} -p ${conda_base}/envs/${conda_env_name} || true; }  # corner case: there is already a conda env with the same name in another conda_base, e.g. NGSAdmin
+  [[ -z "${conda_env_path}" ]] && { mamba env create -n ${conda_recipe_env_name} -f $conda_recipe || true; }
+  [[ -n "${conda_env_path}" && "$arg_force" == true && ! -d "${conda_base}/envs/${conda_env_name}" ]] && { mamba env create -f $conda_recipe -p ${conda_base}/envs/${conda_env_name} || true; }  # corner case: there is already a conda env with the same name in another conda_base, e.g. NGSAdmin
   [[ -n "${conda_env_path}" ]] && logwarn "A Conda environment with the name \"$conda_env_name\" is already present. Skipping environment creation."
   # conda activate $conda_env_name
   set -eu
@@ -271,7 +286,7 @@
   [[ -f $repo_path/reference_db/kraken/hash.k2d ]] && status_kraken="${green}OK${normal}"
   [[ -f $repo_path/reference_db/mash/mashDB.msh ]] && status_mash="${green}OK${normal}"
   [[ -f $repo_path/reference_db/taxonkit/names.dmp ]] && status_taxonkit="${green}OK${normal}"
-  [[ -f ${conda_env_path:-$CONDA_PREFIX}/lib/python3.7/site-packages/quast_libs/busco/busco.py ]] && status_busco="${green}OK${normal}"
+  [[ -d ${conda_env_path:-$CONDA_PREFIX}/lib/python3.7/site-packages/quast_libs/busco/bacteria ]] && status_busco="${green}OK${normal}"
   [[ -f ${conda_env_path:-$CONDA_PREFIX}/lib/python3.7/site-packages/quast_libs/augustus3.2.3/bin/augustus ]] && status_augustus="${green}OK${normal}"
   logheader "Database Status:"
   echo "
@@ -291,8 +306,10 @@
 logentry "Arguments: ${*:-"No arguments provided"}"
 
 ## Workflow
-conda_recipe_env_name=$(head -n1 $repo_path/envs/${conda_recipe} | cut -d' ' -f2)
+conda_recipe_env_name=${conda_recipe_env_name:-$(head -n1 $conda_recipe | cut -d' ' -f2)}
+show_info
 check_conda $conda_recipe_env_name
+[[ "$interactive" == true ]] && interactive_query
 [[ "$arg_mamba" == true ]] && setup_mamba
 [[ "$arg_busco" == true ]] && complete_busco
 [[ "$arg_databases" == true ]] && download_databases

EOF

# Finalize instalation with script
bash $CONDA_PREFIX/opt/aquamis/scripts/aquamis_setup.sh \
  --env_name $(basename $CONDA_PREFIX) \
  --busco --databases --auto
