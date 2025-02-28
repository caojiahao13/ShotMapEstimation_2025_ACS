Data description.

Game contexts:
    - `RealDataAnalysis/Data/jordan/covariates_home.txt`: Game contexts of Michael Jordan.
    - `RealDataAnalysis/Data/curry/covariates.txt`: Game contexts of Stephen Curry.
    - `RealDataAnalysis/Data/james/covariates.txt`: Game contexts of LeBron James.
    - Each file contains the game records of a player with three columns: 
        * "game": Game ID starts from 1; 
        * "home": 1 for home and 2 for away; 
        * "level": 1 for strong opponent and 2 for weak opponent.
Shot locations and outcomes:
    - `RealDataAnalysis/RealDataAnalysis/Data/jordan/newdata.txt`: Shot locations and outcomes of Michael Jordan.
    - `RealDataAnalysis/Data/curry/del_curry.txt`: Shot locations and outcomes of Stephen Curry.
    - `RealDataAnalysis/Data/james/del_james.txt`: Shot locations and outcomes of LeBron James.
    - Each file contains the shot records of a player with columns including:
        * "Game": Game ID starts from 1 which matches "game" values in the game context files; 
        * "Missmade": 0 for missed shot and 1 for made shot; 
        * "Xdist" and "Ydist": $(x,y)$-coordinates of shot locations.