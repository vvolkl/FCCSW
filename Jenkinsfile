pipeline {
    agent any

    stages {
        stage('Build-centos7') {
            agent {
              label "centos7"
            }
            steps {
                sh """
                source /cvmfs/fcc.cern.ch/sw/views/releases/externals/94.2.0/x86_64-centos7-gcc62-opt/setup.sh
                make -j `getconf _NPROCESSORS_ONLN` && make test
                
                """
            }
        }
        stage('Build-slc6') {
            agent {
              label "slc6"
            }
            steps {
                sh """
                source /cvmfs/fcc.cern.ch/sw/views/releases/externals/94.2.0/x86_64-slc6-gcc62-opt/setup.sh
                make -j `getconf _NPROCESSORS_ONLN` && make test
                
                """
            }
          }
    }
}
