
#---then Geant4
    IFS=$'\n'
    for line in `$thisdir/bin/geant4-config --datasets`; do
        unset IFS
        read dataset var value <<< $line;
        export $var=$value
    done
    export G4INSTALL=$(cd $(geant4-config --prefix);pwd)
    export G4EXAMPLES=$G4INSTALL/share/Geant4-$(geant4-config --version)/examples

