Arborescence de l'archive : 
    /simulation :
        - propagation_avec_miroir.py : Propagation avec pompe seule dans un système "Matériau(x) linéaire(s)"- Semiconducteur - Miroir
        - propagation_avec_sonde.py : Expérience pompe-sonde simulée pour un système "Matériau(x) linéaire(s)"- Semiconducteur - Miroir
        - variation_epaisseur_para.py : Expérience pompe-sonde simulée pour plusieurs épaisseurs de semiconducteurs pour un système "Matériau(x) linéaire(s)"- Semiconducteur - Miroir
        - variation_tau.py : Expérience pompe-sonde simulée pour plusieurs temps de recombinaison \tau_c
    /comparaison : 
        - data_experimental/ : Dossier contenant les données expérimentales utilisées pour les comparaisons
            - [XXX]ExperimentalPoint.csv : Pompe-sonde pour [XXX]nm de GaAs
            - BEVEL2_summary.csv : Pompe-sonde réalisée en fonction de l'épaisseur pour un pulse d'amplitude 1.7 V/m
        - amplitude.py : Comparaison entre expérience et simulation de la réflectivité max en fonction de l'amplitude d'entrée
        - epaisseur.py : Comparaison entre expérience et simulation de la réflectivité (max et min) en fonction de l'épaisseur du semiconducteur. 
        - temporel.py : Affichage de la réflectivité en fonction du temps pour plusieurs puissance de pulse (il faut les données). 

Comment exécuter un code ? 
    1. Pour créer les dossiers et télécharger les libraries nécessaires au bon fonctionnement, il faut exécuter setup.sh
    2. Activer le venv qui a été créé
    3. Rester dans le dossier racine puis exécuter le code voulu (ex: python comparaison/epaisseur.py)
    4. Les sorties sont enregistrées dans le dossier results.

Comment modifier des paramètres d'entrée ?
    Les paramètres d'entrée sont indiqués en en-tête des codes de simulation. 
    Il n'y a normalement aucun paramètre à modifier à l'intérieur des fonctions.