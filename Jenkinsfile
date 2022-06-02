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
	    sh 'make hdf5'
	    sh 'make nulib'
	    sh 'make lua'
	    sh 'sed -i 's/eos_table_name.*/eos_table_name=\.\.\/tables\/Hempel_SFHoEOS_rho222_temp180_ye60_version_1\.1_20120817\.h5/' external/NuLib/parameters'
	    dir('external/NuLib'){
	        sh './make_table_example'
		sh 'mv NuLib*.h5 NuLib.h5'
	    }
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
