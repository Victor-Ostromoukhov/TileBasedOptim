## data

Dans data sont situés les fichiers de configuration et le profil MatBuilder, ainsi que les matrices.

## Execution

Pour executer ce programme, après avoir configuré les fichiers il faut faire comme suit :

```
cd optimizer/build
cmake ..
make -j8
mv ./bin/OptimMSE2DPipeLine ~/bin/
cd ScriptCreator/build
cmake ..
make -j8
./bin/ScriptCreator
cd ../../Scripts
./LaunchScript.sh
```
Après avoir compilé l'executable de l'optimiseur, mettre celui dans le bin du home.
Ensuite, compiler le créateur de script et l'executer
Enfin, executer le script
Il faut gaire attention à bien configurer les parametre dans les fichiers de configuration

