# vim: sw=4 ts=4 et
DBIND_ENV=dynamicbind
RELAX_ENV=relax
PKGGREP='accelerate\|attn\|awq\|cuda\|datasets\|deepspeed\|einops\|formers\|gptq\|llava\|numba\|numpy\|nvidia\|optimize\|peft\|pydantic\|scikit\|scipy\|sglang\|torch\|triton\|urllib3\|vllm'


# Set to false => remove pre-existing and recreate the env
reuse_dbind=false
reuse_relax=true

if [ -f ~/bin/conda.sh ]; then
    . ~/bin/conda.sh
fi

# Note: env relax is very old -- I've newer openmm installs that are fairly up-to-date
#       allow relaxation with ani or mace accurate force fields.

if $reuse_dbind && conda env list | grep -q "^${DBIND_ENV} "; then
    echo "GOOD: Reusing env ${DBIND_ENV}"
else
    echo "Removing any old env ${DBIND_ENV}"
    if conda env list | grep -q "^${DBIND_ENV} "; then
        conda env remove -n "${DBIND_ENV}"
    fi
    #rm -rf "${CONDA_PREFIX}/envs/${DBIND_ENV}"    # ??
    echo "Recreating env ${DBIND_ENV}"
    conda create -y -n ${DBIND_ENV} python=3.10
    conda deactivate; conda activate ${DBIND_ENV}

    # cuda 11.7 -> 11.8
    set -x
    conda install -y --override-channels -c conda-forge -c nvidia -c pytorch \
        -c schrodinger -c pyg \
        pytorch torchvision torchaudio pytorch-cuda=11.8 transformers \
        pyg pyyaml biopython rdkit openbabel \
        jupyter ipython jupyterlab ipywidgets nb_conda_kernels \
        openmm openmm-ml openmmtools openmm-torch openmm-plumed openmmforcefields
    pip install pyg_lib torch_scatter torch_sparse torch_cluster torch_spline_conv -f https://data.pyg.org/whl/torch-2.0.0+cu118.html
    pip install e3nn fair-esm spyrmsd
    set +x
fi
conda list -n ${DBIND_ENV} > ${DBIND_ENV}.list
cat ${DBIND_ENV}.list | grep -i "${PKGGREP}"

if $reuse_relax && conda env list | grep -q "^${RELAX_ENV} "; then
    echo "GOOD: Reusing env ${RELAX_ENV}"
else
    echo "Removing any old env ${RELAX_ENV}"
    if conda env list | grep -q "^${RELAX_ENV} "; then
        conda env remove -n "${RELAX_ENV}"
    fi
    #rm -rf "${CONDA_PREFIX}/envs/${RELAX_ENV}"    # ??
    echo "Recreating env ${RELAX_ENV}"
    conda create -y -n ${RELAX_ENV} python=3.8
    conda deactivate; conda activate ${RELAX_ENV}

    conda install -y --override-channels -c conda-forge openmm pdbfixer libstdcxx-ng \
        openmmforcefields openff-toolkit ambertools=22 compilers biopython \
        openmm-torch openmm-ml openmmtools plumed py-plumed # this line added!
fi
conda list -n ${RELAX_ENV} > ${RELAX_ENV}.list
cat ${RELAX_ENV}.list | grep -i "${PKGGREP}"

