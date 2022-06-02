pipeline {
    triggers { pollSCM('') }  // Run tests whenever a new commit is detected.
    agent { dockerfile {args '--gpus all'}} // Use the Dockerfile defined in the root Flash-X directory
    stages {

        //=============================//
    	// Set up submodules and amrex //
        //=============================//
    	stage('Prerequisites'){ steps{
	    sh 'mpicc -v'
	    sh 'nvidia-smi'
	    sh 'nvcc -V'
	    sh 'git submodule update --init'
	    sh 'cp make.inc.template make.inc'
	    sh 'make nulib'
	}}

	stage('0D SR'){ steps{
	    sh 'make clean; DEBUG=1 NDIMS=0 DO_GR=0 make all -j'
	    sh 'exe/MultiDArrayTest'
	    sh 'make -C tests/inelastic_scatter'
	    // sh 'make -C tests/blackbody'
	}}

    } // stages{

    post {
        always {
	    cleanWs(
	        cleanWhenNotBuilt: true,
		deleteDirs: true,
		disableDeferredWipeout: false,
		notFailBuild: true,
		patterns: [[pattern: 'external', type: 'EXCLUDE']] ) // allow amrex to be cached
	}
    }

} // pipeline{
