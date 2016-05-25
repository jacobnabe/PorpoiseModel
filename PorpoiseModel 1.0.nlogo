;  PorpoiseModel version 1.0

;  This is the source code for the individual-based model used for evaluating the
;  cumulative effects of by-catch and noise from wind farms and ships on the harbour porpoise 
;  population in the inner Danish waters. Please refer to the scientific publication for detailed 
;  documentation: Nabe-Nielsen, J., Sibly, R.M., Tougaard, J., Teilmann, J. & Sveegaard, S. (2014) 
;  "Effects of noise and by-catch on a Danish harbour porpoise population". Ecological 
;  Modelling, 272, 242–251.


; The model was created as part of the projects:
; EFFECTS OF WIND FARMS ON HARBOUR PORPOISE BEHAVIOUR AND POPULATION DYNAMICS
; funded by Miljogruppen: By- og Landskabsstyrelsen, Energistyrelsen,
; DONG and Vattenfall A/S
; and 
; BRIDGES AS BARRIERS PORPOISE MODELLING: DOES THE GREAT BELT BRIDGE HINDER 
; MOVEMENT OF HARBOUR PORPOISES IN THE GREAT BELT
; funded by Fehmern Belt A/S
;

; Copyright (C) 2016, Jacob Nabe-Nielsen <jnn@bios.au.dk>
; 
; This program is free software; you can redistribute it and/or modify it 
; under the terms of the GNU General Public License version 2 and only 
; version 2 as published by the Free Software Foundation.
; 
; This program is distributed in the hope that it will be useful,
; but WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
; GNU General Public License for more details.
; 
; You should have received a copy of the GNU General Public License
; along with this program; if not, write to the Free Software
; Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


; The model was developed and tested using NetLogo version 4.1. Development ended 2011-11-29

 
; debug levels: 
;   0   no debugging
;   1   debug porp-avoid-land / water depth problems
;   2   write turning angles before and after avoiding land
;   3   debugging turning angles -- esp. loop control
;   4   debugging attraction vector caused by reference memory
;   5   debugging porp-get-exp-food-val -- expected value of future food
;   6   debugging average energy level in 40x40 km blocks
;   7   debugging porpoise dispersal
;   8   debugging deterrence from ships and wind turbines
;   9   debugging reproduction

; behavioural modes in model: 
;   0   Behaviour produced by simple non-autocorrelated Markov process (except when porps seek to avoid land); 
;       variables calibrated using dead-reckoning data.
;   1   Like model 0, but introducing autocorrelation in turning angles and speed, and sharper turning angles 
;       at low speed (i.e. a correlated random walk (CRW) model). 
;   2   Extends model 1 behaviour by introducing a desire to return to return to areas with food (i.e. a reference
;       memory, cf Van Moorter, Oikos 2009). Food is eaten in the entire cell once the porp has been there. 
;       Food doesn't affect step length (intrinsic behaviour).
;   3   Extends model 2 by introducing long-distance dispersal, energy consumption and a response to wind farms and ships.
;   4   Extends model 3 by introducing mortality and reproduction (and altered energy use when with calf). 
;       Reproduction and mortality occurs daily.


extensions [ gis profiler ]
globals [
  time-step
  sim-day
  month                  ; using 30-day months
  quarter                ; season of the year (used for updating food growth)
  year
  prev-month
  prev-quarter
  prev-year
  path                   ; path to directory for input/output, one dir for each area
  outfile                ; simulation output
  bathy-data 
  block-data             ; block numbers 1-60 (40 x 40 km blocks, numbered by row, starting in upper left corner)
  block-centres-x        ; List of x coordinates for block centres
  block-centres-y        ; List of y coordinates for block centres
  block-values           ; Expected relative quality of blocks, based on MEAN of Maxent values per block. Used for selecting dispersal direction
  list-of-dead-age       ; Age of death for all animals that die. Reset every year
  list-of-dead-day       ; Day of death for all animals that die. Reset every year
  disttocoast-data
  food-prob-data
  maxent-level-data
  mean-maxent-in-quarters ; standardized average maxent level in each quarter
  ; daily-surv-prob          ; Average daily survival in each porpoise age class. From Caswell et al 1998
  food-upd-interval      ; Number of days between updating amount of food. Updated daily in model >=1
  maxU                   ; Maximum utility of a patch, set to 1 here
  wind-turb-list         ; List of wind turbines. Each entry contains a list with id (char), x (num), y (num), and intensity (num)
  aarhus-odden-route     ; Positions on the ship route Aarhus-Odden. List of lists w two elements.
  great-belt-route       ; Positions on the ship T-route through Kattegat. List of lists w two elements.
  kattegat-sound-route   ; Positions on the ship T-route through the Sound. List of lists w two elements.
  age-distrib            ; Number of porps in each age class (ages >=25 y aggregated)
  ship-ao-list           ; List of ships in the Aarhus-Odden route
  ship-gb-list           ; List of ships in the Great-Belt route
  ship-kat-list          ; List of ships in the Kattegat-Sound route
  turb-who-list          ; 
  xllcorner 
  yllcorner
  corr-logmov            ; correlation in movement distance in CRW +
  corr-angle             ; correlation in direction in CRW +
  deploy_x               ; x-coord of porp 0
  deploy_y               ; y-coord of porp 0
  turn-right             ; random variable; turn right if = 1, left if = -1
  min-depth
  movlgt-list            ; list of 100 move lengths for porpoise 0 (for histogram)
  angle-list             ; list of 100 turning angles for porpoise 0 (for histogram)
  max-movlgt             ; monotonously increasing, move length/100m
  mean-movlgt
  memory-max             ; Maximum number of half-hour steps the amount of food can be remembered (120 steps is 2.5 days)
  inertia-const          ; A, set to 0.001 like in ms
  ref-mem-strength-list-fixed   ; replaces ref-mem-strength-list, uses rR = 0.10
  work-mem-strength-list-fixed  ; replaces work-mem-strength-list, uses rW = 0.20
  use-exp-food-val       ; get more attracted to the CRW path if food was found recently
  vt                     ; resultant attraction vector, resulting from reference memory of food availability (model >=2)
  CRW-contrib            ; length of vector pointing in direction predicted by CRW
  ; MR-contrib             ; length of vector pointing in direction of remembered food
 ]
 
breed [ porps porp ]
breed [ turbs turb ]
breed [ ao-ships ship ]
breed [ gb-ships ship ]
breed [ kat-ships ship ]


patches-own [ 
  bathymetry
  block
  disttocoast            ; Distance to coast (in m)
  food-prob              ; randomly distributed patches, 1 cell large (prob = 1 inside patches, 0 outside)
  maxent-level           ; growth rate for the food in the patches (from quarterly MAXENT predictions)
  food-level             ; current amount of food in cell (zero outside food patches)
]

porps-own [ 
  ;lgth                   ; Length in cm
  ;weight                 ; Weight in kg
  age                    ; Age in years (remember, 360 days per year)
  age-of-maturity        ; Age when becoming receptive, in years
  pregnancy-status       ; 0 (unable to mate, young/low energy); 1 (unable to mate, pregnant); 2 (ready to mate)
  mating-day             ; Day of year. Most mating occurs in August (Lockyer 2003)
  ds-mating              ; Days since mating. -99 if not pregnant
  dsg-birth              ; Days since giving birth. -99 if not with lactating calf
  with-lact-calf         ; true/false, with lactating calf
  energy-level           ; Porpoises get energy by eating and loose energy by moving.
  energy-level-sum       ; Sum of energy levels. Reset to 0 every day
  energy-level-daily     ; List with average energy for last ten days. Latest days first.
  disp-type              ; Disperse away from low-energy area. 0 if not dispersing, 1 if dispersing far, 2 if...
  disp-target            ; List with x and y coord of the patch that the porp attempts to disperse to (not UTM)
  prev-angle             ; Last turning angle (not the heading!)
  pres-angle             ; Present turning angle
  prev-logmov            ; Previous Log10 (move length [measured in 100-m steps])
  pres-logmov            ; Present Log10 (move length [measured in 100-m steps])
  enough-water-ahead     ; Turn to avoid land if false
  pos-list               ; Coordinates of previous positions -- one per 30 min, latest positions in left end (= item 0)
  pos-list-daily         ; Coordinates of previous 10 daily positions -- daily, corresponding to energy-level-daily
  ; Vars added in model 2
  ; ref-mem-strength-list  ; Memory decay with time (logistic decrease, function of rR)
  ; work-mem-strength-list ; Memory of poor quality patches, decays with time (logistic decrease, function of rW)
  ; work-mem-updated       ; Whether a list of working memory has been updated
  stored-util-list       ; Up_t ; Remembered feeding success (after memory decay) -- latest positions left
  deter-vt               ; Vector (= list) determining which direction a porp is deterled from wind turbines and ships, and how much
  deter-strength
  deter-time-left
  VE-total               ; total value of food expected to be found in the future
]

ao-ships-own [
  route
  id
  speed
  impact
  past-loc
  next-loc
]

gb-ships-own [
  route
  id
  speed
  impact
  past-loc
  next-loc
]

kat-ships-own [
  route
  id
  speed
  impact
  past-loc
  next-loc
]

turbs-own [
  id                      ; User defined identifier
  impact                  ; Deterrence effect relative to standard Roedsand-turbine
]

; Landscape variables:

to landsc-setup
  ;; (for this model to work with NetLogo's new plotting features,
  ;; __clear-all-and-reset-ticks should be replaced with clear-all at
  ;; the beginning of your setup procedure and reset-ticks at the end
  ;; of the procedure.)
  __clear-all-and-reset-ticks
  reset-timer
  ; Note that setting the coordinate system here is optional, as
  ; long as all of your datasets use the same coordinate system.
  ; gis:load-coordinate-system (word "data/" projection ".prj")
  ; Load datasets:
  set path word "raster-data/" area      ; note that "/" may be a problem on Windows
  set path word path "/"
  set bathy-data gis:load-dataset word path "bathy.asc"
  set food-prob-data gis:load-dataset word path "patches.asc"
  set maxent-level-data gis:load-dataset word path "quarter1.asc"
  set block-data gis:load-dataset word path "blocks.asc"
  set block-centres-x ( list 50 150 250 350 450 550  50 150 250 350 450 550  50 150 250 350 450 550  50 150 250 350 450 550  50 150 250 350 450 550  50 150 250 350 450 550  50 150 250 350 450 550  50 150 250 350 450 550  50 150 250 350 450 550  50 150 250 350 450 550 )
  set block-centres-y ( list 950 950 950 950 950 950 850 850 850 850 850 850 750 750 750 750 750 750 650 650 650 650 650 650 550 550 550 550 550 550 450 450 450 450 450 450 350 350 350 350 350 350 250 250 250 250 250 250 150 150 150 150 150 150  50  50  50  50  50  50 )
  set disttocoast-data gis:load-dataset word path "disttocoast.asc"
  ; Init landscape-related variables
  set food-upd-interval 1
  set maxU 1
  set time-step 0 
  set sim-day 0
  set year 0
  set month 1
  set quarter 1
  
  ; Use Kattegat-area throughout
  set xllcorner 529473
  set yllcorner 5972242

  ; Set the world envelope to the union of all of our dataset's envelopes
  gis:set-world-envelope (gis:envelope-union-of 
    (gis:envelope-of bathy-data)
    (gis:envelope-of food-prob-data)
    (gis:envelope-of maxent-level-data)
    )
  ; This is the preferred way of copying values from a raster dataset
  ; into a patch variable: in one step, using gis:apply-raster.
  gis:apply-raster bathy-data bathymetry
  gis:apply-raster food-prob-data food-prob
  gis:apply-raster maxent-level-data maxent-level
  gis:apply-raster block-data block
  gis:apply-raster disttocoast-data disttocoast
  
  ; set amount of food -- if there is a chance that food is present
  set mean-maxent-in-quarters (list 0.515686364223653 0.888541219760357 0.841346010536882 1) ; setting global variable here 
  ask patches [ ifelse food-prob > 0 [ set food-level maxU * maxent-level / item (quarter - 1) mean-maxent-in-quarters ] [ set food-level 0 ] ]
    
  landsc-display
  print ""
  let tmp-t word "Setup time: " timer
  print word tmp-t " sec"
end


to turbs-import-pos  ; get wind turbine positions and shipping routes
  if (wind-farms != "off") [
    let wf word "wind-farms/" wind-farms
    set wf word wf ".txt"
    let line-txt ""
    let line-lst ( list " " " " " " " ")
    set wind-turb-list ( list line-lst ); a list of line-lst elements
    let line-lgt 0
    if (not file-exists? wf) [ 
      print word "Cannot load file: " wf
      stop
    ]
    file-open wf
    let i 0  ; chars
    let j 0  ; word nbr
    let k 0  ; line nbr
    let one-char ""
    let one-word ""
    let in-word true
    while [ not file-at-end? ] [
      set line-txt file-read-line
      set line-txt word line-txt " "
      set line-lgt length line-txt
      set i 0
      set j 0
      set in-word true
      while [ i < line-lgt and k > 0 ] [  ; first line (k = 0) contains the header, must be omitted
        set one-char substring line-txt i (i + 1)
        if ((one-char = "\n" or one-char = "\t" or one-char = " ") and in-word) [  ; that is, prev. char was part of a word, current isn't
          set in-word false
          if (j > 0) [ set one-word read-from-string one-word ] ; cast to numeric
          set line-lst replace-item j line-lst one-word
          set j j + 1
        ]
        if ((one-char = "\t" or one-char = " ") and (not in-word)) [
          set one-word ""
        ]
        if (one-char != "\t" and one-char != " ") [
          set in-word true
          set one-word word one-word one-char 
        ]
        set i i + 1
      ] 
      if (k > 0) [ set wind-turb-list lput line-lst wind-turb-list ]
      if (k = 1) [ set wind-turb-list remove-item 0 wind-turb-list ]
      set k k + 1
    ]
    file-close
  ; plot the wind turbines
  print word "Showing wind turbines at: " wind-farms
  if (pile-driving) [ print "Pile driving" ]  
  ]
end ; end turbs-import-pos


to turbs-deter-porps
  ; Deterrence by wind farms after construction phase
  let curr-deter 0 ; current amount of deterring
  let deter-dist std-deterrence-dist / 400 ; number of grid-cells where a wind turbine with impact 1 (standard deterrence strength) affects a porpoise
  let turb-pos list -9 -9
  let porp-pos list -9 -9
  let dist-to-turb -9
  if (not pile-driving) [
    ask turbs [
      set turb-pos list ([xcor] of self) ([ycor] of self)
      ask porps in-radius (deter-dist * impact) [  ; i.e. porps <deter-dist away (although porps can hear only turbs <200 away, the dist has to be larger to account for porp jumping)
        set porp-pos list ([xcor] of self) ([ycor] of self)
        set dist-to-turb distancexy (item 0 turb-pos) (item 1 turb-pos)
        set curr-deter ([impact] of myself * deter-dist ) - dist-to-turb       ; deterring-strength decreases linearly with distance to turbine, decreases to 0 at 400 m
        if (deter-strength < curr-deter ) [ ; become deterred if not already more scared of other wind turbine
          set deter-strength curr-deter
          set deter-vt replace-item 0 deter-vt (curr-deter * ((item 0 porp-pos) - (item 0 turb-pos)))   ; vector pointing away from turbine
          set deter-vt replace-item 1 deter-vt (curr-deter * ((item 1 porp-pos) - (item 1 turb-pos)))
          set deter-time-left deter-time ; how long to remain affected
        ]
        if (debug = 8) [ 
          set color 105 
          type word "(porp " who
          type word ") dist-to-turb " myself
          type word ": " round (dist-to-turb * 400)
          print word " m, curr.deter: " round curr-deter
        ]
        ; Porpoises nearby stop dispersing (which could force them to cross over disturbing agents very fast)
        set disp-type 0 
      ]
    ] 
  ]
  ; Deterrence by wind farms during construction
  if (pile-driving) [
    if ((year = 1) or (year = 10) or (year = 20) or (year = 30) or (year = 40)) [
      ; Make sure that only one turbine is constructed at a time
      ; let ddd round (round (sim-day / 2) - (year - 1) * 360) ; two days construction time per turbine, max one year 
      let ddd round (sim-day - (year - 1) * 360) / 2 ; two days construction time per turbine, max one year 
      carefully [
        ask turb (item ddd turb-who-list) [ 
          set color black 
          set turb-pos list ([xcor] of self) ([ycor] of self)
          ask porps in-radius (20000 / 400) [  ; i.e. porps <20 km away get disturbed
            set porp-pos list ([xcor] of self) ([ycor] of self)
            set dist-to-turb distancexy (item 0 turb-pos) (item 1 turb-pos)
            set curr-deter (1 * (20000 / 400) ) - dist-to-turb       ; deterring-strength decreases linearly with distance to turbine, decreases to 0 at 400 m
            if (deter-strength < curr-deter ) [ ; become deterred if not already more scared of other wind turbine
              set deter-strength curr-deter
              set deter-vt replace-item 0 deter-vt (curr-deter * ((item 0 porp-pos) - (item 0 turb-pos)))   ; vector pointing away from turbine
              set deter-vt replace-item 1 deter-vt (curr-deter * ((item 1 porp-pos) - (item 1 turb-pos)))
              set deter-time-left deter-time ; how long to remain affected
                                             ; If very close to construction (< 8 km), move right away
              if (dist-to-turb * 400 < 16000) [
                facexy (item 0 turb-pos) (item 1 turb-pos)
                rt 180
                carefully [ 
                  fd 2 
                  set disp-type 2
                ] [ set disp-type 2 ]
              ]
            ]
            if (debug = 8) [ 
              set color 105 
              type word "(porp " who
              type word ") dist-to-turb " myself
              type word ": " round (dist-to-turb * 400)
              print word " m, curr.deter: " round curr-deter
            ]
            ; Porpoises nearby stop dispersing (which could force them to cross over disturbing agents very fast)
            set disp-type 0            
          ] 
        ]
      ]
      [ print "Runtime error in procedure turbs-deter-porps" ] 
    ]
    
  ]  
end ; end turbs-deter-porps


to turbs-setup
  set turb-who-list [ ]
  ifelse (length wind-turb-list = 0) [ print "No wind turbines plotted" ] [
    let i 0
    create-turbs length wind-turb-list [ 
      ifelse pile-driving [ set color 8 ] [ set color black ]
      set size 3
      set shape "x"
    ]
    ask turbs [
      set id (item 0 (item i wind-turb-list))
      let new-x (( item 1 (item i wind-turb-list) ) - xllcorner ) / 400
      let new-y (( item 2 (item i wind-turb-list) ) - yllcorner ) / 400
      set impact (item 3 (item i wind-turb-list))
      setxy new-x new-y
      set i i + 1
      set turb-who-list lput who turb-who-list
    ]
  ]
end ; end turbs-setup


to file-noise-debug ; setup file and write; Turbines in line have ycor = 706
    let outfile-for-debug "output/Homogeneous/file-noise-debug.txt" ; overwrite existing data
    if (not file-exists? outfile-for-debug) [ 
      file-open outfile-for-debug
      print ( "Writing to output/Homogeneous/file-noise-debug.txt" )
      file-print (" deterr-coef std-deter-dist id time-step utm-x utm-y") ; header line (space-separated). 
      file-close
    ]
    file-open outfile-for-debug
    ask porps [
      ;if ( ycor < (706 + 5) and ycor > (706 - 5)) [
      if ( ycor < (706 + 20) and ycor > (706 - 20)) [
        file-write deterrence-coeff
        file-write std-deterrence-dist
        file-write who
        file-write time-step
        file-write xcor
        file-write ycor
        file-print ""    ; followed by CR
      ]
    ]
    file-close
    
end ; end file-noise-debug



to ships-import  ; get ship impacts and speeds [km / h] in all routes. Separate lists for different routes:
  ; aarhus-odden route: ship-ao-list read from file
  if (incl-ships) [
    let ship-ao "ships/Aarhus-Odden.txt" 
    let ship-gb "ships/Great-Belt.txt" 
    let ship-kat "ships/Kattegat-Sound.txt" 
    let line-txt ""
    let line-lst ( list " "  " "  " " )
    set ship-ao-list ( list line-lst ); a list of line-lst elements
    set ship-gb-list ( list line-lst ); a list of line-lst elements
    set ship-kat-list ( list line-lst ); a list of line-lst elements
    let line-lgt 0
    if (not file-exists? ship-ao) [ 
      print word "Cannot load file: " ship-ao
      stop
    ]
    if (not file-exists? ship-gb) [ 
      print word "Cannot load file: " ship-gb
      stop
    ]
    if (not file-exists? ship-kat) [ 
      print word "Cannot load file: " ship-kat
      stop
    ]
    ; One route at a time -- close file afterwards:
    file-open ship-ao
    let i 0  ; chars
    let j 0  ; word nbr
    let k 0  ; line nbr
    let one-char ""
    let one-word ""
    let in-word true
    while [ not file-at-end? ] [
      set line-txt file-read-line
      set line-txt word line-txt " "
      set line-lgt length line-txt
      set i 0
      set j 0
      set in-word true
      while [ i < line-lgt and k > 0 ] [  ; first line (k = 0) contains the header, must be omitted
        set one-char substring line-txt i (i + 1)
        if ((one-char = "\n" or one-char = "\t" or one-char = " ") and in-word) [  ; that is, prev. char was part of a word, current isn't
          set in-word false
          if (j > 0) [ set one-word read-from-string one-word ] ; cast to numeric
          set line-lst replace-item j line-lst one-word
          set j j + 1
        ]
        if ((one-char = "\t" or one-char = " ") and (not in-word)) [
          set one-word ""
        ]
        if (one-char != "\t" and one-char != " ") [
          set in-word true
          set one-word word one-word one-char 
        ]
        set i i + 1
      ] 
      if (k > 0) [ set ship-ao-list lput line-lst ship-ao-list ]
      set k k + 1
    ]
    file-close
    ; next route:
    file-open ship-gb
    set i 0  ; chars
    set j 0  ; word nbr
    set k 0  ; line nbr
    set one-char ""
    set one-word ""
    set in-word true
    while [ not file-at-end? ] [
      set line-txt file-read-line
      set line-txt word line-txt " "
      set line-lgt length line-txt
      set i 0
      set j 0
      set in-word true
      while [ i < line-lgt and k > 0 ] [  ; first line (k = 0) contains the header, must be omitted
        set one-char substring line-txt i (i + 1)
        if ((one-char = "\n" or one-char = "\t" or one-char = " ") and in-word) [  ; that is, prev. char was part of a word, current isn't
          set in-word false
          if (j > 0) [ set one-word read-from-string one-word ] ; cast to numeric
          set line-lst replace-item j line-lst one-word
          set j j + 1
        ]
        if ((one-char = "\t" or one-char = " ") and (not in-word)) [
          set one-word ""
        ]
        if (one-char != "\t" and one-char != " ") [
          set in-word true
          set one-word word one-word one-char 
        ]
        set i i + 1
      ] 
      if (k > 0) [ set ship-gb-list lput line-lst ship-gb-list ]
      set k k + 1
    ]
    file-close
   ; next route:
    file-open ship-kat
    set i 0  ; chars
    set j 0  ; word nbr
    set k 0  ; line nbr
    set one-char ""
    set one-word ""
    set in-word true
    while [ not file-at-end? ] [
      set line-txt file-read-line
      set line-txt word line-txt " "
      set line-lgt length line-txt
      set i 0
      set j 0
      set in-word true
      while [ i < line-lgt and k > 0 ] [  ; first line (k = 0) contains the header, must be omitted
        set one-char substring line-txt i (i + 1)
        if ((one-char = "\n" or one-char = "\t" or one-char = " ") and in-word) [  ; that is, prev. char was part of a word, current isn't
          set in-word false
          if (j > 0) [ set one-word read-from-string one-word ] ; cast to numeric
          set line-lst replace-item j line-lst one-word
          set j j + 1
        ]
        if ((one-char = "\t" or one-char = " ") and (not in-word)) [
          set one-word ""
        ]
        if (one-char != "\t" and one-char != " ") [
          set in-word true
          set one-word word one-word one-char 
        ]
        set i i + 1
      ] 
      if (k > 0) [ set ship-kat-list lput line-lst ship-kat-list ]
      set k k + 1
    ]
    file-close

  ]
end ; end ships-import


to ships-setup
  let i 0
  if (not (length ship-ao-list < 2)) [ ; aarhus-odden. First line contains header
    set ship-ao-list remove-item 0 ship-ao-list ; remove header
    let n-ao-ships length ship-ao-list
    set aarhus-odden-route (list (list 119 633) (list 157 605) (list 197 597) (list 288 584) )
    create-ao-ships n-ao-ships [
      set route "aarhus-odden"
      set past-loc random length aarhus-odden-route
      ifelse (past-loc = (length aarhus-odden-route - 1)) [set next-loc past-loc - 1] [set next-loc past-loc + 1]
      set shape "arrow"
      set color black
      set size 10
      setxy (item 0 (item past-loc aarhus-odden-route)) (item 1 (item past-loc aarhus-odden-route))
      facexy (item 0 (item next-loc aarhus-odden-route)) (item 1 (item next-loc aarhus-odden-route))
    ]
    ask ao-ships [
      set id item 0 (item i ship-ao-list)
      set speed item 1 (item i ship-ao-list)
      set impact item 2 (item i ship-ao-list)
      set i i + 1
    ]
  ]
  set i 0
  if (not (length ship-gb-list < 2)) [ ; great belt route. First line contains header
    set ship-gb-list remove-item 0 ship-gb-list ; remove header
    let n-gb-ships length ship-gb-list
    set great-belt-route (list (list 288.66 999) (list 366.56 804.34) (list 210.04 539.42) (list 205.78 489.69) (list 249.14 404.9) (list 259.2 369.07) (list 249.49 308.61) (list 225.86 263.82) (list 228.83 227.3) (list 377.59 154.1) (list 438.14 156.74) (list 445.4 168) (list 524.38 261.43) (list 599 315.93))
    create-gb-ships n-gb-ships [
      set route "great-belt"
      set past-loc random length great-belt-route
      ifelse (past-loc = (length great-belt-route - 1)) [set next-loc past-loc - 1] [set next-loc past-loc + 1]
      set shape "arrow"
      set color black
      set size 10
      setxy (item 0 (item past-loc great-belt-route)) (item 1 (item past-loc great-belt-route))
      facexy (item 0 (item next-loc great-belt-route)) (item 1 (item next-loc great-belt-route))
    ]
    ask gb-ships [
      set id item 0 (item i ship-gb-list)
      set speed item 1 (item i ship-gb-list)
      set impact item 2 (item i ship-gb-list)
      set i i + 1
    ]
  ]
  set i 0
  if (not (length ship-kat-list < 2)) [ ; kattegat sound route. First line contains header
    set ship-kat-list remove-item 0 ship-kat-list ; remove header
    let n-kat-ships length ship-kat-list
    set kattegat-sound-route (list (list 288.66 999) (list 366.56 804.34) (list 401.52 678.13) (list 478 628) (list 487.77 620.61) (list 508.29 575) (list 542.26 500.6) (list 506.3 407.47) (list 599 315.93) )
        create-kat-ships n-kat-ships [
      set route "kattegat-sound"
      set past-loc random length kattegat-sound-route
      ifelse (past-loc = (length kattegat-sound-route - 1)) [set next-loc past-loc - 1] [set next-loc past-loc + 1]
      set shape "arrow"
      set color black
      set size 10
      setxy (item 0 (item past-loc kattegat-sound-route)) (item 1 (item past-loc kattegat-sound-route))
      facexy (item 0 (item next-loc kattegat-sound-route)) (item 1 (item next-loc kattegat-sound-route))
    ]
    ask kat-ships [
      set id item 0 (item i ship-kat-list)
      set speed item 1 (item i ship-kat-list)
      set impact item 2 (item i ship-kat-list)
      set i i + 1
    ]
  ]
end

to ships-move
  let the-route "NA"
  let mov-lgt 0

  ask ao-ships [
    set mov-lgt 2.5 * speed / 2  ; km / t to steps / 30 min
    set the-route aarhus-odden-route
    ; approaching next location on route?
    if distancexy (item 0 (item next-loc the-route)) (item 1 (item next-loc the-route)) < mov-lgt [
      ; check wheter position numbers should increase, or if ships should turn around:
      let tmp next-loc
      ifelse (next-loc > past-loc and next-loc < (length(the-route) - 1) ) [ set next-loc next-loc + 1 ] [ set next-loc next-loc - 1 ]
      set past-loc tmp
      if  (next-loc < 0) [ set next-loc 1 ]
      facexy (item 0 (item next-loc the-route)) (item 1 (item next-loc the-route))
    ]
    fd mov-lgt
  ]
  ask gb-ships [
    set mov-lgt 2.5 * speed / 2  ; km / t to steps / 30 min
    set the-route great-belt-route
    ; approaching next location on route?
    if distancexy (item 0 (item next-loc the-route)) (item 1 (item next-loc the-route)) < mov-lgt [
      ; check wheter position numbers should increase, or if ships should turn around:
      let tmp next-loc
      ifelse (next-loc > past-loc and next-loc < (length(the-route) - 1) ) [ set next-loc next-loc + 1 ] [ set next-loc next-loc - 1 ]
      set past-loc tmp
      if  (next-loc < 0) [ set next-loc 1 ]
      facexy (item 0 (item next-loc the-route)) (item 1 (item next-loc the-route))
    ]
    fd mov-lgt
  ]
  ask kat-ships [
    set mov-lgt 2.5 * speed / 2  ; km / t to steps / 30 min
    set the-route kattegat-sound-route
    ; approaching next location on route?
    if distancexy (item 0 (item next-loc the-route)) (item 1 (item next-loc the-route)) < mov-lgt [
      ; check wheter position numbers should increase, or if ships should turn around:
      let tmp next-loc
      ifelse (next-loc > past-loc and next-loc < (length(the-route) - 1) ) [ set next-loc next-loc + 1 ] [ set next-loc next-loc - 1 ]
      set past-loc tmp
      if  (next-loc < 0) [ set next-loc 1 ]
      facexy (item 0 (item next-loc the-route)) (item 1 (item next-loc the-route))
    ]
    fd mov-lgt
  ]
end


to ships-deter-porps
  let curr-deter 0 ; current amount of deterring
  let deter-dist std-deterrence-dist / 400 ; number of grid-cells where a wind turbine or ship with impact 1 (standard deterrence strength) affects a porpoise
  let ship-pos list -9 -9
  let porp-pos list -9 -9
  let dist-to-ship -9
  ask ao-ships [
    set ship-pos list ([xcor] of self) ([ycor] of self)
    ask porps in-radius (deter-dist * impact) [  ; i.e. porps <deter-dist away (although porps can hear only ships <200 away, the dist has to be larger to account for porp jumping)
      set porp-pos list ([xcor] of self) ([ycor] of self)
      set dist-to-ship distancexy (item 0 ship-pos) (item 1 ship-pos)
      set curr-deter ([impact] of myself * deter-dist ) - dist-to-ship       ; deterring-strength decreases linearly with distance to turbine, decreases to 0 at 400 m
      if (deter-strength < curr-deter ) [ ; become deterred if not already more scared of other wind turbine
        set deter-strength curr-deter
        set deter-vt replace-item 0 deter-vt (curr-deter * ((item 0 porp-pos) - (item 0 ship-pos)))   ; vector pointing away from turbine
        set deter-vt replace-item 1 deter-vt (curr-deter * ((item 1 porp-pos) - (item 1 ship-pos)))
        set deter-time-left deter-time ; how long to remain affected
      ]
      if (debug = 8) [ 
        set color 105 
        type word "who: " who
        type word " dist-to-ship " myself
        print word ": " dist-to-ship 
      ]
      ; Porpoises nearby stop dispersing (which could force them to cross over disturbing agents very fast)
      set disp-type 0 
    ] 
  ]  
  ; next ship breed
  set curr-deter 0 ; current amount of deterring
  set deter-dist std-deterrence-dist / 400 ; number of grid-cells where a wind turbine or ship with impact 1 (standard deterrence strength) affects a porpoise
  set ship-pos list -9 -9
  set porp-pos list -9 -9
  set dist-to-ship -9
  ask gb-ships [
    set ship-pos list ([xcor] of self) ([ycor] of self)
    ask porps in-radius (deter-dist * impact) [  ; i.e. porps <deter-dist away (although porps can hear only ships <200 away, the dist has to be larger to account for porp jumping)
      set porp-pos list ([xcor] of self) ([ycor] of self)
      set dist-to-ship distancexy (item 0 ship-pos) (item 1 ship-pos)
      set curr-deter ([impact] of myself * deter-dist ) - dist-to-ship       ; deterring-strength decreases linearly with distance to turbine, decreases to 0 at 400 m
      if (deter-strength < curr-deter ) [ ; become deterred if not already more scared of other wind turbine
        set deter-strength curr-deter
        set deter-vt replace-item 0 deter-vt (curr-deter * ((item 0 porp-pos) - (item 0 ship-pos)))   ; vector pointing away from turbine
        set deter-vt replace-item 1 deter-vt (curr-deter * ((item 1 porp-pos) - (item 1 ship-pos)))
        set deter-time-left deter-time ; how long to remain affected
      ]
      if (debug = 8) [ 
        set color 105 
        type word "who: " who
        type word " dist-to-ship " myself
        print word ": " dist-to-ship 
      ]
      ; Porpoises nearby stop dispersing (which could force them to cross over disturbing agents very fast)
      set disp-type 0 
    ] 
  ]  
  ; next ship breed
  set curr-deter 0 ; current amount of deterring
  set deter-dist std-deterrence-dist / 400 ; number of grid-cells where a wind turbine or ship with impact 1 (standard deterrence strength) affects a porpoise
  set ship-pos list -9 -9
  set porp-pos list -9 -9
  set dist-to-ship -9
  ask kat-ships [
    set ship-pos list ([xcor] of self) ([ycor] of self)
    ask porps in-radius (deter-dist * impact) [  ; i.e. porps <deter-dist away (although porps can hear only ships <200 away, the dist has to be larger to account for porp jumping)
      set porp-pos list ([xcor] of self) ([ycor] of self)
      set dist-to-ship distancexy (item 0 ship-pos) (item 1 ship-pos)
      set curr-deter ([impact] of myself * deter-dist ) - dist-to-ship       ; deterring-strength decreases linearly with distance to turbine, decreases to 0 at 400 m
      if (deter-strength < curr-deter ) [ ; become deterred if not already more scared of other wind turbine
        set deter-strength curr-deter
        set deter-vt replace-item 0 deter-vt (curr-deter * ((item 0 porp-pos) - (item 0 ship-pos)))   ; vector pointing away from turbine
        set deter-vt replace-item 1 deter-vt (curr-deter * ((item 1 porp-pos) - (item 1 ship-pos)))
        set deter-time-left deter-time ; how long to remain affected
      ]
      if (debug = 8) [ 
        set color 105 
        type word "who: " who
        type word " dist-to-ship " myself
        print word ": " dist-to-ship 
      ]
      ; Porpoises nearby stop dispersing (which could force them to cross over disturbing agents very fast)
      set disp-type 0 
    ] 
  ]  
end ; end ships-deter-porps


to update-block-values
  ; Update relative values for blocks 1-60 every quarter, depending on location -- based on MEAN of maxent values for block
  let block-val-homo ( list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 )
  let block-val-kat-1 ( list 0 0.007 0.0042 0.0022 8e-04 0 0 0.0454 0.0522 0.0176 0.0185 0 0 0.0458 0.0348 0.0062 0.0182 0.4446 0 0.2898 0.1013 0.0769 0.126 0.4615 0.4439 0.5104 0.5741 0.4593 0.6447 0.6498 0.8685 0.5063 0.971 0.1697 0.27 0.1261 0.9438 0.7346 1 0.3994 0.1576 0.2355 0.7718 0.9082 0.8524 0.2161 0.4916 0.2975 0.8397 0.7533 0.7199 0.885 0.727 0.1057 0 0 0.8284 0.7116 0.4866 0 )
  let block-val-kat-2 ( list 0 0.5443 0.4624 0.5383 0.6492 0 0 0.5012 0.4789 0.4938 0.6435 0 0 0.5225 0.5076 0.4741 0.5791 0.5663 0 0.9995 0.6703 0.6761 0.8228 0.6815 0.871 0.9514 0.8679 0.5161 0.6802 0.7584 0.9344 0.9618 1 0.3191 0.2393 0.2157 0.815 0.813 0.9247 0.5042 0.212 0.1856 0.769 0.7131 0.711 0.2889 0.1534 0.118 0.7809 0.5221 0.5176 0.374 0.244 0.1908 0 0 0.4939 0.4498 0.5051 0 )
  let block-val-kat-3 ( list 0 0.5465 0.4406 0.4862 0.4404 0 0 0.4853 0.3784 0.4794 0.5802 0 0 0.567 0.4035 0.4414 0.5272 0.6176 0 1 0.5381 0.382 0.6403 0.6813 0.8465 0.9411 0.8096 0.5441 0.7883 0.7735 0.9179 0.8931 0.9637 0.3427 0.3232 0.3117 0.9205 0.8392 0.8785 0.4564 0.3102 0.305 0.8802 0.7424 0.6202 0.322 0.2576 0.2078 0.9401 0.5889 0.5387 0.4164 0.2382 0.1777 0 0 0.5698 0.4275 0.4118 0 )
  let block-val-kat-4 ( list 0 0.4315 0.3406 0.2776 0.5347 0 0 0.5089 0.3525 0.259 0.3855 0 0 0.6569 0.4462 0.3708 0.4415 0.7523 0 0.8659 0.5877 0.4377 0.6658 0.8918 0.8696 0.8685 0.9511 0.8375 0.7861 0.8468 0.9827 0.9134 0.932 0.6272 0.6809 0.5617 0.9968 1 0.9454 0.8803 0.6762 0.4845 0.9085 0.9043 0.8077 0.6977 0.5416 0.4581 0.8676 0.7344 0.7909 0.7304 0.64 0.6211 0 0 0.8921 0.8097 0.8367 0 )
  if ( area = "Homogeneous" ) [ set block-values block-val-homo ]
  if ( area = "Kattegat" and quarter = 1 ) [ set block-values block-val-kat-1 ]
  if ( area = "Kattegat" and quarter = 2 ) [ set block-values block-val-kat-2 ]
  if ( area = "Kattegat" and quarter = 3 ) [ set block-values block-val-kat-3 ]
  if ( area = "Kattegat" and quarter = 4 ) [ set block-values block-val-kat-4 ]
end

  
to landsc-display  ; Updates the display variable
  no-display
  if (disp-var = "bathymetry") [ 
    let min-bathymetry gis:minimum-of bathy-data
    let max-bathymetry gis:maximum-of bathy-data
    if (max-bathymetry - min-bathymetry = 0) [ set max-bathymetry max-bathymetry + 0.01 ] ; to avoid err in homogeneous landsc
    ask patches
    [ ; note the use of the "<= 0 or >= 0" technique to filter out 
      ; "not a number" values, as discussed in the documentation.
      set pcolor 39
      if (bathymetry <= 0) or (bathymetry >= 0)
      [ set pcolor scale-color blue bathymetry max-bathymetry (min-bathymetry + 3) ] 
    ]
  ]
  
  if (disp-var = "maxent-level") [ ; potential rate of food increase, from quarterly MAXENT preditions
    let min-food-gr gis:minimum-of maxent-level-data
    let max-food-gr gis:maximum-of maxent-level-data
    if (max-food-gr - min-food-gr = 0) [ set max-food-gr max-food-gr + 0.01 ] ; to avoid err in homogeneous landsc
    ask patches
    [ 
      set pcolor 39
      if (maxent-level <= 0) or (maxent-level >= 0)
;      [ set pcolor scale-color green maxent-level  ((log max-food-gr 10) + 1) ((log (min-food-gr + 0.1) 10) + 0.6)  ] 
      [ set pcolor scale-color green maxent-level  1 -0.2 ] ; should work for probabilities
    ]
  ]

  if (disp-var = "food-prob") [ ; randomly distributed food patches
    let min-food-prob gis:minimum-of food-prob-data
    let max-food-prob gis:maximum-of food-prob-data
    ask patches
    [ 
      set pcolor 39
      if (food-prob <= 0) or (food-prob >= 0)
      [ set pcolor scale-color violet food-prob max-food-prob min-food-prob ] 
    ]
  ]

  if (disp-var = "food-level") [ ; actual amount of food in patches
    let min-food-level 0
    let max-food-level maxU
    ask patches
    [ 
      set pcolor 39
      if (maxent-level <= 0) or (maxent-level >= 0) [
        ;[ set pcolor scale-color green food-level (1.0 * max-food-level) (1.0 * min-food-level) ] 
        if ( food-level = 0 ) [ set pcolor white ]
        if ( food-level > 0 and food-level <= 0.02 * maxU ) [ set pcolor 45 ]
        if ( food-level > 0.02 and food-level <= 0.1 * maxU ) [ set pcolor 27 ]
        if ( food-level > 0.1 * maxU and food-level <= 0.25 * maxU ) [ set pcolor 66 ]
        if ( food-level > 0.25 * maxU and food-level <= 0.5 * maxU ) [ set pcolor 63 ]
        if ( food-level > 0.5 * maxU and food-level < 0.99 * maxU ) [ set pcolor 61 ]
        if ( food-level > 0.99 * maxU ) [ set pcolor 61 ]
      ]
    ]
  ]
  if ( disp-var = "blocks" ) [ 
    let min-block 1
    let max-block 60
    ask patches
    [
      set pcolor 39
      if (block <= 0) or (block >= 0)
      [ set pcolor scale-color red block max-block min-block ] 
    ]
  ]
  display
end  ; end landsc-display


to landsc-upd-food
  ; maxent-level is patch specific, between 0 and 1 (MAXENT-based); food-growth-rate (rU) is global variable
  ask patches with [ food-prob > 0 and food-level < (maxU * maxent-level) ] [  
    let f-lev food-level + ( food-growth-rate * food-level * ( 1 - food-level / (maxU * maxent-level / item (quarter - 1) mean-maxent-in-quarters) ) )
    if (abs (f-lev - food-level) > 0.001) [
      repeat 47 [ set f-lev f-lev + ( food-growth-rate * f-lev * ( 1 - f-lev / (maxU * maxent-level / item (quarter - 1) mean-maxent-in-quarters) ) ) ]
    ]  ; If the food level is really low, let food grow 48 times -- like growing every half-hour step, only faster
    set food-level f-lev
    ; here maxent-level is MAXENT prediction and food-growth-rate is a universal calibrated variable
  ]
end

to landsc-upd-maxent-level-map
  if ( area = "Kattegat" ) [
    let yyy word "quarter" quarter
    set yyy word yyy ".asc"
    set maxent-level-data gis:load-dataset word path yyy
    gis:apply-raster maxent-level-data maxent-level
    landsc-display
  ]
end


; Porpoise variables and methods

to porps-setup-ref
  ; reference porpoises (with who = 0) -- deployment information
  if (area = "Kattegat") [
;    set deploy_x ( 619399.6 - xllcorner ) / 400  ; start-position
;    set deploy_y ( 6358048 - yllcorner ) / 400
    ; E of Sealand
    set deploy_x 527  
    set deploy_y 397
    ; N Kattegat:
    set deploy_x 218  
    set deploy_y 955
  ]
  if (area = "Homogeneous") [
    set deploy_x ( 619399.6 - xllcorner ) / 400  ; start-position, in pixels
    set deploy_y ( 6148048 - yllcorner ) / 400
  ]
  setxy deploy_x deploy_y
  set color red
  print "Porps deployed (red dot = porp 0)"
end

to porps-setup
  set memory-max 120
  set inertia-const 0.001
  set corr-logmov 0.94
  set corr-angle 0.26
  set vt list 0 0
  set ref-mem-strength-list-fixed ( list 0.999 0.9989 0.9988 0.9987 0.9985  0.9984  0.9982  0.9981  0.9979  0.9976  0.9974  0.9972  0.9969  0.9966  0.9962  0.9958  0.9954  0.995  0.9945  0.9939  0.9933  0.9926  0.9919  0.9911  0.9902  0.9893  0.9882  0.987  0.9858  0.9843  0.9828  0.9811  0.9793  0.9772  0.975  0.9726  0.9699  0.967  0.9638  0.9603  0.9565  0.9523  0.9478  0.9428  0.9375  0.9316  0.9252  0.9183  0.9108  0.9027  0.8939  0.8844  0.8742  0.8632  0.8514  0.8387  0.8252  0.8108  0.7954  0.7792  0.7619  0.7438  0.7248  0.7048  0.684  0.6624  0.64  0.617  0.5934  0.5692  0.5447  0.5199  0.4949  0.4699  0.445  0.4203  0.396  0.3721  0.3487  0.326  0.304  0.2828  0.2626  0.2432  0.2248  0.2074  0.1909  0.1755  0.161  0.1475  0.1349  0.1233  0.1125  0.1025  0.0933  0.0848  0.0771  0.0699  0.0634  0.0575  0.0521  0.0471  0.0426  0.0386  0.0349  0.0315  0.0284  0.0257  0.0232  0.0209  0.0189  0.017  0.0153  0.0138  0.0125  0.0112  0.0101  0.0091  0.0082  0.0074 )
  set work-mem-strength-list-fixed ( list 0.9990  0.9988  0.9986  0.9983  0.9979  0.9975  0.9970  0.9964  0.9957  0.9949  0.9938  0.9926  0.9911  0.9894  0.9873  0.9848  0.9818  0.9782  0.9739  0.9689  0.9628  0.9557  0.9472  0.9372  0.9254  0.9116  0.8955  0.8768  0.8552  0.8304  0.8022  0.7705  0.7351  0.6962  0.6539  0.6086  0.5610  0.5117  0.4617  0.4120  0.3636  0.3173  0.2740  0.2342  0.1983  0.1665  0.1388  0.1149  0.0945  0.0774  0.0631  0.0513  0.0416  0.0336  0.0271  0.0218  0.0176  0.0141  0.0113  0.0091  0.0073  0.0058  0.0047  0.0037  0.0030  0.0024  0.0019  0.0015  0.0012  0.0010  0.0008  0.0006  0.0005  0.0004  0.0003  0.0003  0.0002  0.0002  0.0001  0.0001  0.0001  0.0001  0.0001  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 )
  let age-dist-list (list 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  4  4  4  4  4  4  4  4  4  4  4  4  4  4  4  4  4  4  4  4  4  4  5  5  5  5  5  5  5  5  5  5  5  5  5  5  6  6  6  6  6  6  6  6  6  6  6  6  6  6  6  6  6  6  7  8  8  8  8  8  8  8  8  9  9  9  9  9  9  9  9  9  9  9  9 10 10 10 10 10 10 11 11 11 11 11 12 12 12 12 12 12 12 13 13 13 13 14 14 14 14 15 15 15 15 18 18 19 19 21 22) ;stranded + bycaught animals in Lockyer & Kinze 2003 (NAMMCO)
  ;let birth-rate (list 0 0 0.136 0.417 0.818 0.714 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833)
  set max-movlgt 0
  set mean-movlgt 0
  set use-exp-food-val false
  set CRW-contrib -9999
  ; set MR-contrib -9999 
  set min-depth 1               ; water depth where porps very rarely swim (according to Jakob)
  set turn-right 1
  set movlgt-list list (0) (0)  ; two inputs required...  list used for histogramming
  set movlgt-list remove-item 0 movlgt-list
  set movlgt-list remove-item 0 movlgt-list
  set angle-list list (0) (0)  ; two inputs required...  list used for histogramming
  set list-of-dead-age [ ]
  set list-of-dead-day [ ] 
  set age-distrib [ ]
  set angle-list remove-item 0 angle-list
  set angle-list remove-item 0 angle-list
  create-porps n-porps
  ask porps [ 
      set age item (random (length age-dist-list)) age-dist-list      
      set mating-day 0  ; set in yearly-tasks
      set ds-mating -99
      set dsg-birth -99
      ifelse ( age > age-of-maturity )[ set pregnancy-status 2 ][ set pregnancy-status 0 ]
      if ( pregnancy-status = 2 and random-float 1 < 0.68 ) [  ;  become pregnanat with prob. taken from Read & Hohn 1995
          set pregnancy-status 1
          set ds-mating 360 - round( random-normal (7.5 * 360 / 12) 20 )
      ]
      set with-lact-calf false
      set age-of-maturity 3.44  ;see Read (1990) and Caswell for age of maturity.
      set energy-level random-normal 10 1
      set energy-level-sum 0
      set energy-level-daily ( list 0 0 0 0 0 0 0 0 0 0 )
      set disp-type 0
      set disp-target list 0 0
      set prev-logmov 0.8 ; unit = patchsize, i.e. times 100 m
      set prev-angle 10
      let jj 0
      let in-water false
      while [ jj < 100 and in-water = false ] [
        set deploy_x random-xcor
        set deploy_y random-ycor
        setxy deploy_x deploy_y
        if ( [ bathymetry ] of patch-here > 1 ) [ set in-water true ]
        set jj jj + 1
      ]
      set enough-water-ahead true
      let pos list deploy_x deploy_y
      set pos-list ( list pos )
      let tmp (list 0 0)
      set pos-list-daily ( list tmp tmp tmp tmp tmp tmp tmp tmp tmp tmp )
      ; set ref-mem-strength-list [ ]
      ; set work-mem-strength-list [ ]
      ; set work-mem-updated false
      set VE-total 0
      set deter-vt list 0 0
      set deter-strength 0
      set deter-time-left 0
      set stored-util-list  [ ]
      set color orange
      set shape "circle"
      set size 4
      if (who = 0) [ porps-setup-ref ]
      set pen-size 0.05
  ]
  carefully [ 
    set movlgt-list lput ( 10 ^ ( [ prev-logmov ] of porp 0 ) ) movlgt-list
  ]
  [
    set movlgt-list lput 1 movlgt-list
  ]
  if ( not (is-turtle? ( porp 0 )) ) [ porps-setup ]    ; strange that I have to do this... the porpoise isn't always created in the first go
  
  landsc-display ; update displayed amount of food etc
end

to porps-upd-deter
  ask porps [
    if (deter-time-left = 0) [
      set deter-strength 0
      set deter-vt replace-item 0 deter-vt 0
      set deter-vt replace-item 1 deter-vt 0
    ]
    if (deter-time-left > 0) [
      ; reduce deterrence by half every time step
      set deter-time-left deter-time-left - 1
      set deter-strength deter-strength / 2  ; 
      set deter-vt replace-item 0 deter-vt ( item 0 deter-vt / 2 )
      set deter-vt replace-item 1 deter-vt ( item 1 deter-vt / 2 )
    ]
  ]
end ; end porps-upd-deter

to porp-check-depth
  ; Check that there is enough water at all steplengths ahead, set enough-water-ahead to false if < min-depth
  set enough-water-ahead true
  let pres-mov ( 10 ^ pres-logmov )                                                                           ; because pres-logmov may have changed in the porp-avoid-land procedure
  let dd ceiling ( pres-mov / 0.1 )                                                                           ; number of 10-m steps to check water depth at
  let ee 0
  let depth-list list (0) ( [ bathymetry ] of patch-ahead pres-mov )
  set depth-list remove-item 0 depth-list
  repeat dd [
    set ee ee + 1
    set depth-list lput ( [ bathymetry ] of patch-ahead ( ee * 0.1 ) ) depth-list
  ]
  if ( not (( length depth-list ) = length ( filter [ ? > 0 ] depth-list )) ) [                               ; i.e. if some depths on the list aren't > 0
    set enough-water-ahead false
  ]
end


to porp-avoid-land
  ; If shallow water ahead, turn right or left depending on where water is deeper. Turn as little as possible.
  ; Don't do the turning here, but change angle to be turned in porp-std-move or porp-markov-mov.
  ; Note that the emergency procedure "avoid-beh 5" is found in porp-std-move
  let rand-ang random 10
  let avoid-beh 0
  let pres-mov ( 10 ^ pres-logmov )
  let bath-l [ bathymetry ] of patch-left-and-ahead (40 + rand-ang) pres-mov
  let bath-r [ bathymetry ] of patch-right-and-ahead (40 + rand-ang) pres-mov
  ; alternative kinds of evasive behaviour: 
  if ( bath-r >= min-depth or bath-l >= min-depth ) [
    set avoid-beh 1  ; evasive behaviour type 1
    ifelse ( bath-r >= min-depth and bath-l >= min-depth ) 
      [ ifelse ( bath-r >= bath-l ) ; comparison can be true only if neither bath-r or bath-l are NaN, i.e. if both are > min-depth
        [ set pres-angle pres-angle + (40 + rand-ang) ]
        [ set pres-angle pres-angle - (40 + rand-ang) ]
      ]
      [ ifelse ( bath-r >= min-depth ) 
        [ set pres-angle pres-angle + (40 + rand-ang) ]
        [ set pres-angle pres-angle - (40 + rand-ang) ]
      ]
  ]
  ; else try turning more aprubtly ( = 70 deg )
  if not ( bath-r >= min-depth or bath-l >= min-depth ) [
    set avoid-beh 2  ; evasive behaviour type 2
    set bath-l [ bathymetry ] of patch-left-and-ahead (70 + rand-ang) pres-mov
    set bath-r [ bathymetry ] of patch-right-and-ahead (70 + rand-ang) pres-mov
    if ( bath-r >= min-depth or bath-l >= min-depth ) [
      ifelse ( bath-r >= min-depth and bath-l >= min-depth ) 
        [ ifelse ( bath-r >= bath-l ) ; comparison can be true only if neither bath-r or bath-l are NaN, i.e. if both are > min-depth
          [ set pres-angle pres-angle + (70 + rand-ang) ]
          [ set pres-angle pres-angle - (70 + rand-ang) ]
        ]
        [ ifelse ( bath-r >= min-depth ) 
          [ set pres-angle pres-angle + (70 + rand-ang) ]
          [ set pres-angle pres-angle - (70 + rand-ang) ]
        ]
    ]
  ]
  if not ( bath-r >= min-depth or bath-l >= min-depth ) [
    set avoid-beh 3  ; evasive behaviour type 3
    set bath-l [ bathymetry ] of patch-left-and-ahead (120 + rand-ang) pres-mov
    set bath-r [ bathymetry ] of patch-right-and-ahead (120 + rand-ang) pres-mov
    if ( bath-r >= min-depth or bath-l >= min-depth ) [
      ifelse ( bath-r >= min-depth and bath-l >= min-depth ) 
        [ ifelse ( bath-r >= bath-l ) ; comparison can be true only if neither bath-r or bath-l are NaN, i.e. if both are > min-depth
          [ set pres-angle pres-angle + (120 + rand-ang) ]
          [ set pres-angle pres-angle - (120 + rand-ang) ]
        ]
        [ ifelse ( bath-r >= min-depth ) 
          [ set pres-angle pres-angle + (120 + rand-ang) ]
          [ set pres-angle pres-angle - (120 + rand-ang) ]
        ]
    ]
  ]  
  if not ( bath-r >= min-depth or bath-l >= min-depth ) [
    ; if everything else fails, turn around
    set avoid-beh 4  ; evasive behaviour type 4
    let j 0
    porp-check-depth
    while [ not enough-water-ahead and j < length pos-list ] [
      facexy (item 0 (item j pos-list)) (item 1 (item j pos-list))  ; each item on pos-list contains a list with a x and a y-coordinate
      setxy (item 0 (item j pos-list)) (item 1 (item j pos-list))
      set j j + 1
      porp-check-depth
      if (j = 20) [ set enough-water-ahead true ]
    ]    
  ]
  if ( debug = 1 ) [ 
    let tmp-list list ("beh =") avoid-beh 
    set tmp-list lput ("; tck =") tmp-list
    set tmp-list lput time-step tmp-list
    write tmp-list 
    let tmp7 word "; " round pres-angle
    print word tmp7 " degr." 
  ]
end ; end porp-avoid-land


to porp-markov-move
  ; Movements based on dead-reckoning data -- first calc distance, then turning angle
  set pres-logmov ( 0.5 + random-normal 0 0.25 ) 
  let pres-mov ( 10 ^ pres-logmov )
  set pres-angle random-normal 0 40
  if ( abs pres-angle > 60 ) [ set pres-angle (1 + random-float 0.5) * pres-angle ]  ; make angle dist more leptokurtic
  right pres-angle
  ;
  ; Turn to avoid swimming on land if necessary:
  ; ( section copied in from porp-std-move)
  let dd ceiling ( pres-mov / 0.25 )  ; number of 25-m steps to check water depth at
  let goto-avoid-land false
  if (not ( [ bathymetry ] of patch-ahead pres-mov >= min-depth ) ) [ set goto-avoid-land true ]
  repeat dd [
    if ( not ( [ bathymetry ] of patch-ahead ( dd * 0.25 ) >= min-depth ) ) [ set goto-avoid-land true ]  ; must use "not >= " rather than " < " for catching NaN's
    set dd dd - 1
  ]
  if ( goto-avoid-land ) [ porp-avoid-land ]
  set pres-mov ( 10 ^ pres-logmov )  ; because pres-logmov may have changed in the porp-avoid-land procedure
;  let ready-to-move true
  ; test again:
  set dd ceiling ( pres-mov / 0.1 )  ; number of 10-m steps to check water depth at
  let ee 0
  let depth-list list (0) ( [bathymetry] of patch-ahead pres-mov )
  set depth-list remove-item 0 depth-list
  repeat dd [
    set ee ee + 1
    set depth-list lput ( [ bathymetry ] of patch-ahead ( ee * 0.1 ) ) depth-list
  ]
  if ( not (( length depth-list ) = length ( filter [ ? > 0 ] depth-list )) ) [ ; i.e. if some items on the list aren't < 0
    uphill bathymetry
    if ( debug = 1 ) [ 
      show word "Tick = " time-step
      show word "Moved to deeper patch, depth = " ([bathymetry] of patch-here) 
    ]
  ]
  ;
  ; move
;  if (ready-to-move) [ 
  fd pres-mov / 4  ; divide by four when cell size is 400 x 400 m
;     ]
  ; Remember current moves for the next iteration
  set pres-logmov log pres-mov 10
  set prev-angle pres-angle
  set prev-logmov pres-logmov
end


to porp-ref-mem-turn
  ; Move towards places visited previously if food was found there and they aren't too far away or forgotten.
  let bb ( [food-level] of patch-here )  ; Stationary food species. The stored intrisic patch utility for t=0. Initially it is either 0, 1, 
  ; or -9999, but grows logistically after food is eaten
  set vt list -9999 -9999  ; init to catch errors
  if not ( abs(bb) >= 0 ) [  
    ; There are errors in food availability -- sometimes Na is calculated even though depth is > 0. Catch error here
    set bb 0
    if ( debug = 4 ) [ 
      print "Replaced NaN food value with 0"
      print patch-here
    ]
  ]
  set stored-util-list fput bb stored-util-list
  
  ; Update reference memory strength for past locations
  ; Replaced by "ref-mem-strength-list-fixed"
  let max-mem 0.999
  ; set ref-mem-strength-list fput max-mem ref-mem-strength-list   ; Nearly perfect memory of what just happened, cf Van Moorter et al 2009
  ; let ii 1
  ; while [ ii < length ref-mem-strength-list ] [  
  ;   let MRPt item (ii - 1) ref-mem-strength-list  ; (reference memory for patch p at time t)
  ;   let reduced-mem MRPt - ref-mem-decay * (1 - MRPt) * MRPt ; logistic decrease in memory with time, see eqn 2 in Van Moorter
  ;   set ref-mem-strength-list replace-item ii ref-mem-strength-list reduced-mem
  ;   set ii ii + 1
  ; ]

  ; Set patch value for each past location -- perceived patch utility (= reference memory x intrinsic patch utility (stuff eaten)), divided by DIRECT distance
  let perceived-util-list [ ]
  set perceived-util-list lput (item 0 stored-util-list * max-mem) perceived-util-list
  let tmp list (0) (0)
  let attr-vector-list [ ] ; each element in the list is a list with an x-and a y-direction. Vector for first element (this place) has length 0, the others have length 1
  set attr-vector-list lput tmp attr-vector-list
  let one-attr-vector [ ]
  let vector-lgt 0
  let ii 1
  let dist-to-foodpos 0
  while [ ii < length pos-list ] [
    if (item ii stored-util-list = 0) [  ; save time by skipping dist measure when there is no food -- changed 101214
      set perceived-util-list lput 0 perceived-util-list
      set one-attr-vector list 0 0
      set attr-vector-list lput one-attr-vector attr-vector-list
    ] 
    if (not ( item ii stored-util-list = 0 )) [  ; save time by skipping dist measure when there is no food -- changed 101214
      set dist-to-foodpos (distancexy ( item 0 (item ii pos-list) ) ( item 1 (item ii pos-list) ))
      ifelse (dist-to-foodpos < 1E-20 )
        [ set perceived-util-list lput 9999 perceived-util-list ]      ; arbitrary large value for close dist
        [ set perceived-util-list lput ( (item ii stored-util-list) * (item ii ref-mem-strength-list-fixed) / dist-to-foodpos ) perceived-util-list ]
      ; = utility * memory / distance
      ; Create attraction vectors; unit-vectors pointing towards the patches in memory
      set one-attr-vector list ((item 0 (item ii pos-list)) - xcor)  ((item 1 (item ii pos-list)) - ycor)
      ; make sure that it works with wrapping landscapes:
      if ( item 0 one-attr-vector > world-width / 2 ) [ set one-attr-vector replace-item 0 one-attr-vector ( item 0 one-attr-vector - world-width ) ]
      if ( item 0 one-attr-vector < (- world-width / 2 ) ) [ set one-attr-vector replace-item 0 one-attr-vector ( item 0 one-attr-vector + world-width ) ]
      if ( item 1 one-attr-vector > world-height / 2 ) [ set one-attr-vector replace-item 1 one-attr-vector ( item 1 one-attr-vector - world-height ) ]
      if ( item 1 one-attr-vector < (- world-height / 2 ) ) [ set one-attr-vector replace-item 1 one-attr-vector ( item 1 one-attr-vector + world-height ) ]
      set vector-lgt  sqrt( item 0 one-attr-vector * item 0 one-attr-vector + item 1 one-attr-vector * item 1 one-attr-vector )
      if vector-lgt = 0 [
        if ( debug = 4 ) [ 
          show word "attr-vector-lgt = " vector-lgt
          print "skipping to next porp"
        ]
        stop
      ]
      set one-attr-vector replace-item 0 one-attr-vector ((item 0 one-attr-vector) / vector-lgt)
      set one-attr-vector replace-item 1 one-attr-vector ((item 1 one-attr-vector) / vector-lgt)
      set attr-vector-list lput one-attr-vector attr-vector-list
    ]
    set ii ii + 1
  ]
  ; Calculate resultant attraction vector vt as sum of products of individual values and attraction vectors (eqn 5). May have length != 1
  set ii 1  ; no attraction to current pos (t=0)
  let vt-x 0
  let vt-y 0
  while [ ii < length pos-list ] [
    set vt-x vt-x + item ii perceived-util-list * item 0 ( item ii attr-vector-list )
    set vt-y vt-y + item ii perceived-util-list * item 1 ( item ii attr-vector-list )
    set ii ii + 1
  ]
  if ( debug = 4 ) [ 
    type word "Food here: " bb
    type ",  Attr.v: "
    let attr-vect list vt-x (",")
    set attr-vect lput vt-y attr-vect
    print attr-vect
    if (not ( abs(vt-x) >= 0)) [  ; catch missing values
      write "Perc.util: "
      print perceived-util-list
    ]
  ]
  set vt list vt-x vt-y
  
  ; Remove items in distant past to increase execution speed
  ; if ( length ref-mem-strength-list-fixed > memory-max ) [ set ref-mem-strength-list remove-item memory-max ref-mem-strength-list ]
  if ( length stored-util-list > memory-max ) [ set stored-util-list remove-item memory-max stored-util-list ]
end  ; end porp-ref-mem-turn


to porp-work-mem-turn
   ; Influences direction moved in std-move through vector 'vt'
   ; This procedure MUST be called after porp-ref-mem-turn, as it uses the stored-util-list (experienced food at patch) calculated there,
   ; and because it updates the vt vector which is later used in std.move (adds to the previously calculated vt).

  ; Update working memory strength (short-term memory) for past locations
  let max-mem 0.999
  let MWPt 0
  ; set work-mem-strength-list fput max-mem work-mem-strength-list   ; Nearly perfect memory of what just happened, cf Van Moorter et al 2009


;  print work-mem-strength-list ; ### TMP
  let ii 1
;  if ( work-mem-updated = false ) [
;    while [ ii < length work-mem-strength-list ] [  
;      set MWPt item (ii - 1) work-mem-strength-list  ; (working memory for patch p at time t)
;      let reduced-mem MWPt - work-mem-decay * (1 - MWPt) * MWPt ; logistic decrease in memory with time, see eqn 2 in Van Moorter 
;      set work-mem-strength-list replace-item ii work-mem-strength-list reduced-mem
;      set ii ii + 1
;    ]
;  ]
;  set work-mem-updated true
  ; Presently no need to multiply with stored-util-list to get perceived utility -- the animal knows that all food is eaten in the patch.
;  print work-mem-strength-list ; ### TMP


  let tmp list (0) (0)
  let deter-vector-list [ ] ; each element in the list is a list with an x-and a y-direction. Vector for first element (0; this place) has length 0, the others have length 1
  set deter-vector-list lput tmp deter-vector-list 
;  show word "deter-vector-list: " deter-vector-list
  let one-deter-vector [ ]
  let vector-lgt 0
  set ii 1
;  let dist-to-foodpos 0
  while [ ii < length pos-list ] [
    ; Create deterrence vectors; unit-vectors pointing towards the patches in memory
    set one-deter-vector list ((item 0 (item ii pos-list)) - xcor)  ((item 1 (item ii pos-list)) - ycor)
    ; make sure that it works with wrapping landscapes:
    if ( item 0 one-deter-vector > world-width / 2 ) [ set one-deter-vector replace-item 0 one-deter-vector ( item 0 one-deter-vector - world-width ) ]
    if ( item 0 one-deter-vector < (- world-width / 2 ) ) [ set one-deter-vector replace-item 0 one-deter-vector ( item 0 one-deter-vector + world-width ) ]
    if ( item 1 one-deter-vector > world-height / 2 ) [ set one-deter-vector replace-item 1 one-deter-vector ( item 1 one-deter-vector - world-height ) ]
    if ( item 1 one-deter-vector < (- world-height / 2 ) ) [ set one-deter-vector replace-item 1 one-deter-vector ( item 1 one-deter-vector + world-height ) ]
    set vector-lgt sqrt ( item 0 one-deter-vector * item 0 one-deter-vector + item 1 one-deter-vector * item 1 one-deter-vector )
    if vector-lgt = 0 [
      if ( debug = 5 ) [ 
        show word "deter-vector-lgt = " vector-lgt
        print "Haven't moved, skipping to next porp"
      ]
      stop
    ]
    set one-deter-vector replace-item 0 one-deter-vector ((item 0 one-deter-vector) / vector-lgt)
    set one-deter-vector replace-item 1 one-deter-vector ((item 1 one-deter-vector) / vector-lgt)
    set deter-vector-list lput one-deter-vector deter-vector-list
    set ii ii + 1
  ]

  ; Calculate resultant deterrence vector vtd as sum of products of individual values and deterrence vectors
  set ii 1  ; no deterrence from current pos (t=0)
  let vtd-x 0
  let vtd-y 0
  while [ ii < length pos-list ] [
    set vtd-x vtd-x + inertia-const * item ii work-mem-strength-list-fixed * item 0 ( item ii deter-vector-list )
    set vtd-y vtd-y + inertia-const * item ii work-mem-strength-list-fixed * item 1 ( item ii deter-vector-list )
    set ii ii + 1
  ]

  if ( debug = 5 ) [ 
;    print word "work-mem: " work-mem-strength-list
    type "Deter.v: "
    let deter-vect list vtd-x (",")
    set deter-vect lput vtd-y deter-vect
    print deter-vect
    if (length pos-list > 1) [ print word "pos. before: " item 1 pos-list ]
    print word "pos. now: " item 0 pos-list
    print ""
    ; Checked -- works, at least with length pos-list = 2
  ]
  
  ; vtd points towards the previous position, must be subtracted from vt
  set vt replace-item 0 vt ( item 0 vt - vtd-x )
  set vt replace-item 1 vt ( item 1 vt - vtd-y )

; Remove items in distant past to increase execution speed
;  if ( length work-mem-strength-list > memory-max ) [ set work-mem-strength-list remove-item memory-max work-mem-strength-list ]
end ; end porp-work-mem-turn


to porp-get-exp-food-val
  ; Calculate the expaected value (VE-total) of the food to be found in the future based on food found in recent positions x the working memory
  ; Uses the values of the patches in "stored-util-list", calculated in porp-ref-mem-turn

  ; Update working memory strength (short-term memory) for past locations
  let max-mem 0.999
  let MWPt 0
;  set work-mem-strength-list fput max-mem work-mem-strength-list   ; Nearly perfect memory of what just happened, cf Van Moorter et al 2009
;  let ii 1
;  if ( work-mem-updated = false ) [  ; list may have been updated in porp-work-mem-turn
;    while [ ii < length work-mem-strength-list ] [ 
;      set MWPt item (ii - 1) work-mem-strength-list  ; (working memory for patch p at time t)
;      let reduced-mem MWPt - work-mem-decay * (1 - MWPt) * MWPt ; logistic decrease in memory with time, see eqn 2 in Van Moorter 
;      set work-mem-strength-list replace-item ii work-mem-strength-list reduced-mem
;      set ii ii + 1
;    ]
;  ]
;  set work-mem-updated true

  let ii 1
  set VE-total 0
  let max-i min list ( length work-mem-strength-list-fixed ) ( length stored-util-list )
  while [ ii < max-i ] [ 
    set VE-total VE-total + item (ii - 1) work-mem-strength-list-fixed * item (ii - 1) stored-util-list
    set ii ii + 1
  ]

  if ( debug = 5 ) [ 
    print ""
    ; print word "stored-util-list: " stored-util-list
    ; print word "work-mem-strength-list: " work-mem-strength-list
    let ttt (list "VE-total for porp " who ":")
    show word ttt VE-total
  ]
end ; end porp-get-exp-food-val


to porp-std-move
  ; Movements based on dead-reckoning data:
  ; global vars: corr-logmov 0.94 and corr-angle 0.26
  
  ; ### turning angle
  let prev-mov 10 ^ prev-logmov
  let pres-heading heading
  set pres-angle 999
  let j 1
  let tmp-angle 0
  ifelse ( prev-angle < 0 ) [ set tmp-angle prev-angle - 24 ] [ set tmp-angle prev-angle + 24 ]  ; for increasing mean turning angle
  while [ abs (pres-angle) > 180 ]  [ 
    set pres-angle ( tmp-angle * (- corr-angle) + random-normal 0 38 )                  ; Autoreg can't be used for estimating param. as estimated turns are changed if on shallow water. 
    set j j + 1
    if (j = 200) [
      set pres-angle ( pres-angle * 90 / (abs pres-angle))
      if ( debug = 3 ) [ 
        print word "exiting loop 1, ang=" pres-angle
      ]
    ]
  ]
  let sign 0
  ifelse pres-angle < 0 [ set sign -1 ] [ set sign 1 ]
  set pres-angle abs pres-angle ; add the sign again later
  ; Make angle decrease linearly with mov-dist
  let go-on true
  set j 1
  let rnd 0
  while [ go-on ]  [ 
    set rnd random-normal 96 28      ; draws the number to be added to pres-angle
    if ( prev-mov <= 5.50 ) [ 
      set pres-angle pres-angle + rnd - ( rnd * prev-mov / 5.50 )
    ]
    if ( pres-angle < 180 ) [ set go-on false ]  ; remember that turning angle is unsigned here
    set j j + 1
    if (j = 200) [
      set pres-angle ( random 20 + 90 )
      set go-on false
      if ( debug = 3 ) [ 
        print word "exiting loop 2, ang=" pres-angle
      ]
    ]
  ]
  ; if ( abs pres-angle > 55 and abs pres-angle < 180 ) [ set pres-angle (1 - random-float 0.32) * pres-angle ]  ; make turning angle dist more leptokurtic
  set pres-angle pres-angle * sign
  let angle-before-avoid-land pres-angle ; for printing later using debug 2
  right pres-angle
  let angle-turned-right pres-angle ; for updating prev-angle at end of porp-std-move
  set pres-angle 0

  ; ### distance
  set pres-logmov 999
  let porp-max-dist 1.18                                                                                        ; log10 ( max distance a porpoise can travel per half-hour )
  set j 1
  while [ pres-logmov > porp-max-dist ] [ 
    set pres-logmov ( corr-logmov * prev-logmov + random-normal 0.42 0.48 ) 
    set j j + 1
    if (j = 200) [
      if (pres-angle = 0) [set pres-angle pres-angle + 0.00001]
      set pres-angle ( pres-angle * 90 / (abs pres-angle))
      if ( debug = 3 ) [ 
        print word "exiting loop 3, ang=" pres-angle
      ]

    ]    
  ]  
  let pres-mov ( 10 ^ pres-logmov )                                                                             ; This is what is plotted in the histogram
  ;
  ; Turn to avoid swimming on land if necessary:
  set enough-water-ahead false
  let count-i 0
  while [ not enough-water-ahead ] [
    porp-check-depth
    if (not enough-water-ahead) [ porp-avoid-land ]
    set pres-mov ( 10 ^ pres-logmov )                                                      ; because pres-logmov may have changed in the porp-avoid-land procedure
    right pres-angle                                                                       ; angle to turn -- pres-angle -- is changed in porp-avoid-land
    set angle-turned-right (angle-turned-right + pres-angle)
    if (angle-turned-right > 180) [ set angle-turned-right angle-turned-right - 360 ]
    if (angle-turned-right < -180) [ set angle-turned-right angle-turned-right + 360 ]
    set pres-angle 0
    set count-i count-i + 1
    if (count-i = 100) [ 
      set enough-water-ahead true 
      if ( debug = 1 ) [ 
        print "caught water-ahead loop"
      ]
    ]
  ]
  ; test depth again, avoid-beh = 5:
  porp-check-depth
  if (not enough-water-ahead) [
    let prev-heading heading
    let p max-one-of neighbors [ bathymetry ]
    carefully [ face p ] [ ]
    set angle-turned-right (angle-turned-right + pres-angle)
    if (angle-turned-right > 180) [ set angle-turned-right angle-turned-right - 360 ]
    if (angle-turned-right < -180) [ set angle-turned-right angle-turned-right + 360 ]
    set pres-mov 1                                                                                            ; move 100 m towards deeper patch
    if ( debug = 1 ) [ 
      let tmp5 list "beh =  5 ; tck " time-step
      write tmp5 
    ]
  ]
  ;
  ; slow down if turning sharply:
  if ( pres-mov > 10 and (abs angle-turned-right) > 90 ) [ set pres-mov  pres-mov / 5  ] 
  if ( pres-mov > 7 and (abs angle-turned-right) > 50 ) [ set pres-mov  pres-mov / 2  ] 
  ;
  ; Change direction if attracted / deterred by certain areas (model >= 2)
  let total-dx 0
  let total-dy 0
  if ( not use-exp-food-val ) [
    set total-dx (dx * pres-mov) + (item 0 vt)                ; vt isn't used in porp-std-move till here
    set total-dy (dy * pres-mov) + (item 1 vt)                ; note that dx is change in x if taking ONE step forward
    facexy (xcor + total-dx) (ycor + total-dy)
  ]
  if ( use-exp-food-val and model < 3  ) [
    set CRW-contrib inertia-const + pres-mov * VE-total       ; length of vector pointing in direction predicted by CRW (VE-total and pres-mov are porp variables)
    ; set MR-contrib sqrt ( (item 0 vt) * (item 0 vt) + (item 1 vt) * (item 1 vt) )     ; length of vector pointing in direction of remembered food
    set total-dx (dx * CRW-contrib) + (item 0 vt)
    set total-dy (dy * CRW-contrib) + (item 1 vt)
    facexy (xcor + total-dx) (ycor + total-dy)                ; really not needed, it already points that way
  ]
  ; deterrence behaviour -- get scared away from ships and wind turbines
  if ( use-exp-food-val and model >= 3 ) [                    
    set CRW-contrib inertia-const + pres-mov * VE-total       
    set total-dx (dx * CRW-contrib) + (item 0 vt) + ((item 0 deter-vt) * deterrence-coeff)
    set total-dy (dy * CRW-contrib) + (item 1 vt) + ((item 1 deter-vt) * deterrence-coeff)
    facexy (xcor + total-dx) (ycor + total-dy)             
  ]
  
  ; Store turn for calc of turning angle in next step:
  ; let total-turn heading - pres-heading   ; total change in heading, including all adjustments till here. 'pres-heading' was calc in beginning of porp-std-move
  let total-turn subtract-headings heading pres-heading   ; total change in heading, including all adjustments till here. 'pres-heading' was calc in beginning of porp-std-move
  
  ;
  ; Move: 
  ; In the population model all movement lengths are still calculated in 100 m steps, but the cell size has increased from 100 m to 400 m
  ; The step should therefore be divided by 4
  fd pres-mov / 4  ; movement length isn't affected by presence of food
  ;
  if ( debug = 2 ) [ 
    if ( time-step = 0 ) [ print "dist angle-before-avoid-land angle-turned-right x y" ]
    let tmp-var2 (round ( (10 ^ prev-logmov) * 100) / 100)    ; THIS IS IMPORTANT -- the porp turns before it moves, so turning angle is affected by previous moving dist
    set tmp-var2 word tmp-var2 " "
    set tmp-var2 word tmp-var2 angle-before-avoid-land
    set tmp-var2 word tmp-var2 " "
    set tmp-var2 word tmp-var2 angle-turned-right
    set tmp-var2 word tmp-var2 " "
    set tmp-var2 word tmp-var2 ( xcor * 100 + xllcorner )
    set tmp-var2 word tmp-var2 " "
    set tmp-var2 word tmp-var2 ( ycor * 100 + yllcorner )
    print tmp-var2
  ]
  if ( debug = 5 ) [ 
    print word "CRW-contrib: " list ( (dx * (inertia-const + pres-mov * VE-total)) ) ( (dy * (inertia-const + pres-mov * VE-total)) )
    print word "MR-contrib: " vt
    print word "dx, dy (after): " list (total-dx) (total-dy)
    print word "heading (after): " heading
    print word "total-turn: " heading
  ]   
  ;
  ; Remember current moves for the next iteration
  ; if attraction to food alters the movement angle (i.e. vt != 0), this isn't remembered for next step
  ; set prev-angle angle-turned-right  ; so the additional turn due to attraction to food does not influence turning angle in next step 
  set prev-angle total-turn  ; so the additional turn due to attraction to food DOES influence turning angle in next step 
  set prev-logmov log pres-mov 10  ; total steplength, resulting from vt + pres-mov
  ;
  ; test depth one last time, avoid-beh = 6 - move back on same track:
  if (not ([ bathymetry ] of patch-here > 0) ) [
    let prev-heading heading
    if (length pos-list > 1 ) [ facexy (item 0 item 1 pos-list) (item 1 item 1 pos-list) ] ; debugging 110323
    set angle-turned-right (angle-turned-right + pres-angle)
    if (angle-turned-right > 180) [ set angle-turned-right angle-turned-right - 360 ]
    if (angle-turned-right < -180) [ set angle-turned-right angle-turned-right + 360 ]
    setxy (item 0 item 1 pos-list) (item 1 item 1 pos-list)                              ; move 100 m towards deeper patch
    if ( debug = 1 ) [ 
      ; print "; -> "
      let tmp6 list "beh =  6 � ; tck " time-step
      set tmp6 word tmp6 " ; "
      set tmp6 word tmp6 angle-turned-right
      print word tmp6 " degr."
    ]
  ]
  ; update position list:
  let pres-pos list xcor ycor
  set pos-list fput pres-pos pos-list
  if ( length pos-list > memory-max ) [ set pos-list remove-item memory-max pos-list ]   
end ;  end porp-std-move


to porp-upd-energetic-status
  ; 1. Reduce food in the patch that the porp just has left. The amount eaten decreases linearly as the porp's energy level increases from 10 to 20 (=max e)
  ;    this does not affect the porpoise's perception of the quality of the area, and therefore the movement is unaffected.
  ; 2. Adjust porpoise energy level based on amount of food found and time spent per half hour
  ; Increase food level in cells with food-level > 0 AFTERWARDS in order to calc. stored-util-list correctly.
  let food-eaten 0
  let fract-of-food-to-eat 0
  if (energy-level < 20) [ set fract-of-food-to-eat (( 20 - energy-level ) / 10) ]  
  if (fract-of-food-to-eat > 0.99) [ set fract-of-food-to-eat 0.99 ]

  ; remove food from patches:
  ask patch (item 0 (item 1 pos-list)) (item 1 (item 1 pos-list)) [ ; item 0 pos-list was the last added element, i.e. the current position
    if food-level > 0 [ 
      set food-eaten food-level * fract-of-food-to-eat
      set food-level food-level - food-eaten
      if ( food-level < 0.01 ) [ set food-level 0.01 ]                                    ; The minimum food level has a strong impact on how fast food gets back
      if ( food-level > 0 and food-level <= 0.02 * maxU ) [ set pcolor 45 ]
      if ( food-level > 0.02 and food-level <= 0.1 * maxU ) [ set pcolor 27 ]
      if ( food-level > 0.1 * maxU and food-level <= 0.25 * maxU ) [ set pcolor 66 ]
      if ( food-level > 0.25 * maxU and food-level <= 0.5 * maxU ) [ set pcolor 63 ]
      if ( food-level > 0.5 * maxU and food-level < 0.99 * maxU ) [ set pcolor 61 ]
      if ( food-level = 0.99 * maxU ) [ set pcolor 61 ]
    ]
  ] 
  set energy-level energy-level + food-eaten

  ; Scale e-use depending on season and lactation
  let scaling-factor 1
  ; Animals have approx 30% lower energy consumption when the water is cold, Nov-Mar, and approx. 15% lower energy consumption in Oct+Apr (Lockyer et al 2003. Monitoring growth and energy utilization of the harbour porpoise (Phocoena phocoena) in human care. Harbour porpoises in the North Atlantic 5:143-175.)
  if (month = 4 or month = 10) [ set scaling-factor 1.15 ]
  if (month > 4 and month < 10) [ set scaling-factor 1.30 ]
  ; Food consumption increases approx 40% when lactating, there is apparently no effect of pregnancy. (Magnus Wahlberg <magnus@fjord-baelt.dk>, unpubl. data)
  if (with-lact-calf) [ set scaling-factor 1.4 * scaling-factor ]

  ; Prob of dying increases with decreasing energy level
  let yearly-surv-prob (1 - (m-mort-prob-const * exp(- energy-level * x-survival-prob-const) ))
  let step-surv-prob 0
  if (energy-level > 0) [ set step-surv-prob exp( ln(yearly-surv-prob) / (360 * 48) )]  ; called every half hour
  if ( random-float 1 > step-surv-prob ) [ 
    if (not with-lact-calf or energy-level <= 0) [
      set list-of-dead-age lput (floor age) list-of-dead-age
      set list-of-dead-day lput (floor sim-day) list-of-dead-day
      die 
    ]
    ; Better abandoning calf than dying
    if (with-lact-calf) [
      set with-lact-calf false
    ]
  ]
  set energy-level energy-level - ( 0.001 * scaling-factor * e-use-per-30min + (( 10 ^ prev-logmov ) * 0.001 * scaling-factor * e-use-per-km / 0.4 ) )
  set energy-level-sum energy-level-sum + energy-level

  ; Un-weaned calves have an elevated risk of dying that depends on the energy level of their parents (not good at finding food, cannot dive deep, higher predation risk)
  if (with-lact-calf) [
    let yearly-juv-surv-prob yearly-surv-prob ; (1 - (m-juv-mort-const * exp(- energy-level * x-survival-prob-const) ))   ; Assume same k as for adults
    let daily-juv-surv-prob 0
    if (energy-level > 0) [ set daily-juv-surv-prob exp( ln(yearly-juv-surv-prob) / 360 )]

  ]
end ; end porp-upd-energetic-status
 

to porp-move
  if ( model = 0 ) [ porp-markov-move ]
  if ( model = 1 ) [ porp-std-move ]
  if ( model = 2 ) [ 
    set use-exp-food-val true
    porp-ref-mem-turn       ; get attracted to places where food was found. Influences direction moved in std-move through vector 'vt' (global var)
    porp-get-exp-food-val   ; determines the tendency to move following CRW behaviour based on foraging success in recent past
    porp-std-move           ; this is where the porp moves forward
    porp-upd-energetic-status           ; food level increases in 'go' -- affect the landscape and energetic status of the porpoise
  ]
  if ( model >= 3 ) [ 
    set use-exp-food-val true
    porp-ref-mem-turn       ; get attracted to places where food was found. Influences direction moved in std-move through vector 'vt' (global var)
    porp-get-exp-food-val   ; determines the tendency to move following CRW behaviour based on foraging success in recent past
    porp-std-move           ; this is where the porp moves forward and responds to noise by turning away
    porp-upd-energetic-status  ; transform food to energy and spend energy based on step length. Food level in patches increases in 'go'
                            ; mortality and pregnancy status is set in daily-tasks for models > 4
  ]
  if not ( [ bathymetry ] of patch-here > 0 ) [ 
    follow-me
    beep
    user-message "Error, no water"
  ]
end ; end porp-move


to porp-disp-target-select  
  ; deciding where to disperse to based on knowledge of other blocks (each block is 100 x 100 cells = 40 000 x 40 000 m)
  let nbr 0
  let block-quality ( list 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 )
;  if (debug = 7 and (who = 0 or who = 1)) [ print block-quality ]
  while [nbr < length (block-quality) ] [
    let one-dx ( abs (xcor - item nbr block-centres-x) )   ; x dist from porp to block centres -- distancexy doesn't work here
    let one-dy ( abs (ycor - item nbr block-centres-y) ) 
    let block-dist round ( sqrt( one-dx * one-dx + one-dy * one-dy ) )
;    if (debug = 7 and (who = 0 or who = 1)) [ 
;      write (nbr + 1)
;      print word "block-dist: " block-dist 
;    ]
    set block-quality replace-item nbr block-quality ( item nbr block-values )
    ; set quality to 0 for blocks that are too close or too far (divide by 0.4 to convert from km to cells)
    if (block-dist < min-dist-to-target / 0.4 ) [ set block-quality replace-item nbr block-quality 0 ]  
    set nbr nbr + 1
  ]
;  if (debug = 7 and (who = 0 or who = 1)) [ print block-quality ]

  ; numbers of the blocks with highest quality (blocks numbered 0-59, best block first)
  let hi-q-blocks [ ]
  if (n-disp-targets = 4) [ set hi-q-blocks (list 0 0 0 0 ) ]
  if (n-disp-targets = 6) [ set hi-q-blocks (list 0 0 0 0 0 0 ) ]
  if (n-disp-targets = 8) [ set hi-q-blocks (list 0 0 0 0 0 0 0 0) ]
  if (n-disp-targets = 10) [ set hi-q-blocks (list 0 0 0 0 0 0 0 0) ]
  if (n-disp-targets = 12) [ set hi-q-blocks (list 0 0 0 0 0 0 0 0 0 0) ]
  let block-quality-srt sort block-quality   ; quality of the blocks. Last element corresponds to first element in hi-q-blocks
  let i 0
  while [i < length(hi-q-blocks)] [
    let hi-qual-i item ((length block-quality-srt) - i - 1) block-quality-srt  ; block w highest value
    let j 0
    while [ j < length(block-values) ] [
      if (item j block-quality = hi-qual-i) [ 
        set hi-q-blocks replace-item i hi-q-blocks j  ; number of best blocks in left end
      ]
      set j j + 1
    ]
    set i i + 1
  ]
  
  ; select block at random from the twelve blocks with highest quality (where qual = mean.food / dist)
  ; Make sure that porps far north do not try to disperse west
  let the-nbr -9
  set the-nbr random (length hi-q-blocks)
  let sel-block item the-nbr hi-q-blocks
  if ( [ block ] of patch-here < 20 and ( sel-block = 30 or sel-block = 31 or sel-block = 36 or sel-block = 37 or sel-block = 43) ) [ set sel-block sel-block + 2 ]  ; find block slightly more to the east if currently far north
  set disp-target list ( item sel-block block-centres-x ) ( item sel-block block-centres-y )  ; coordinates of cell to move towards
  
  if ( debug = 7 and (who = 0 or who = 1) ) [ 
    print " "
    let tmp word "Disp from: " [ block ] of patch-here
    set tmp word tmp " -> " 
    show word tmp sel-block
  ]
 end ; end porp-disp-target-select


to porp-disp-1  ; disperse to new region
  ; When the porp get closer to land than min-dist-to-land it tries to turn away, else it shifts to disp-2
  let min-dist-to-land 2000 

  ; If porp is SE of Sealand, don't do directed dispersal
  if (xcor > 438 and ycor < 478 ) [
  ;if ([block] of patch-here = 35 or [block] of patch-here = 41 ) [
    set disp-type 2 
    stop
  ] 
  
  ; If porp is N of Djursland or in Little Belt, don't do directed dispersal
  if ([block] of patch-here = 14 or [block] of patch-here = 31 ) [
    set disp-type 2 
    stop
  ] 

  ; Go north if north of Fyn / Funen
  if ([block] of patch-here = 32 ) [
    set heading 0
    set disp-type 2 
    fd 1
    stop
  ] 

  if ( debug = 7 and who = 0 ) [ 
      ask ( porp 0 ) [ pen-down ]
  ]  
  if ( debug = 7 and who = 1 ) [ 
      ask ( porp 1 ) [ pen-down ]
  ]  
  ; Find distance to target block centre
  let the-dx (item 0 disp-target) - xcor 
  let the-dy (item 1 disp-target) - ycor
  let the-dist round ( sqrt( the-dx * the-dx + the-dy * the-dy ) )
  facexy ( xcor + the-dx / 2 ) ( ycor + the-dy / 2 )
;  if ( debug = 7 and (who = 0 or who = 1) ) [ print word "Heading A: " heading ]
      
  ; adjust angle to swim towards deep areas
  let bathymetry-ahead ( list
    -999 ;[ bathymetry ] of patch-at-heading-and-distance (heading - 40) (mean-disp-dist)
    -999 ;[ bathymetry ] of patch-at-heading-and-distance (heading - 30) (mean-disp-dist)
    [ bathymetry ] of patch-at-heading-and-distance (heading - 20) (mean-disp-dist)
    [ bathymetry ] of patch-at-heading-and-distance (heading - 10) (mean-disp-dist)
    [ bathymetry ] of patch-at-heading-and-distance (heading) (mean-disp-dist)
    [ bathymetry ] of patch-at-heading-and-distance (heading + 10) (mean-disp-dist)
    [ bathymetry ] of patch-at-heading-and-distance (heading + 20) (mean-disp-dist)
    -999 ;[ bathymetry ] of patch-at-heading-and-distance (heading + 30) (mean-disp-dist)
    -999 ;[ bathymetry ] of patch-at-heading-and-distance (heading + 40) (mean-disp-dist)
  )
  
  let bathymetry-far-ahead ( list
    -999 ;[ bathymetry ] of patch-at-heading-and-distance (heading - 40) (mean-disp-dist * 8)
    [ bathymetry ] of patch-at-heading-and-distance (heading - 30) (mean-disp-dist * 8)
    [ bathymetry ] of patch-at-heading-and-distance (heading - 20) (mean-disp-dist * 8)
    [ bathymetry ] of patch-at-heading-and-distance (heading - 10) (mean-disp-dist * 8)
    [ bathymetry ] of patch-at-heading-and-distance (heading) (mean-disp-dist * 8)
    [ bathymetry ] of patch-at-heading-and-distance (heading + 10) (mean-disp-dist * 8)
    [ bathymetry ] of patch-at-heading-and-distance (heading + 20) (mean-disp-dist * 8)
    [ bathymetry ] of patch-at-heading-and-distance (heading + 30) (mean-disp-dist * 8)
    -999 ;[ bathymetry ] of patch-at-heading-and-distance (heading + 40) (mean-disp-dist * 8)
  )

  ; Turn up to 20 degr towards deepest water, provided that there is no land further away in that direction
  let i 0
  let good-heading ( list 0 0 0 0 0 0 0 0 0 )
  while [ i < length good-heading ] [
    ifelse ( item i bathymetry-far-ahead > 0 ) [ set good-heading replace-item i good-heading true ] [ set good-heading replace-item i good-heading false ]
    set i i + 1
  ]

  let angles (list -40 -30 -20 -10 0 10 20 30 40 )
  let bathymetry-choice -999
  set i 0
  let sel-angle 0

  set bathymetry-choice max bathymetry-ahead
  while [ i < length bathymetry-ahead ] [ 
    if ( item i bathymetry-ahead = bathymetry-choice ) [ 
      set sel-angle item i angles]
    set i i + 1
  ]
  
  set heading heading + sel-angle
 
;  if (debug = 7 and (who = 0 or who = 1) ) [ 
;    type word "Heading B: " heading 
;    print word ", angle: " sel-angle 
;  ]

  ; Turn to areas far from land if there is land ahead
  let disttocoast-ahead ( list
    -999 ; [ disttocoast ] of patch-at-heading-and-distance (heading - 40) (mean-disp-dist * 2)
    [ disttocoast ] of patch-at-heading-and-distance (heading - 30) (mean-disp-dist * 2)
    [ disttocoast ] of patch-at-heading-and-distance (heading - 20) (mean-disp-dist * 2)
    [ disttocoast ] of patch-at-heading-and-distance (heading - 10) (mean-disp-dist * 2)
    [ disttocoast ] of patch-at-heading-and-distance (heading) (mean-disp-dist * 2)
    [ disttocoast ] of patch-at-heading-and-distance (heading + 10) (mean-disp-dist * 2)
    [ disttocoast ] of patch-at-heading-and-distance (heading + 20) (mean-disp-dist * 2)
    [ disttocoast ] of patch-at-heading-and-distance (heading + 30) (mean-disp-dist * 2)
    -999 ; [ disttocoast ] of patch-at-heading-and-distance (heading + 40) (mean-disp-dist * 2)
  )
  
  ; make sure that there is also water far away
  set i 0
  while [ i < length(disttocoast-ahead) ] [ 
    if ( not item i good-heading ) [ set disttocoast-ahead replace-item i disttocoast-ahead -999 ]
    set i i + 1 
  ]
  
  let disttocoast-choice -999
  let low-water-ahead not (item 1 bathymetry-ahead > min-disp-depth and item 2 bathymetry-ahead > min-disp-depth and item 3 bathymetry-ahead > min-disp-depth and item 4 bathymetry-ahead > min-disp-depth and item 5 bathymetry-ahead > min-disp-depth) ; that is, if some pos ahead are on land
  if ( low-water-ahead or [ disttocoast ] of patch-here < min-dist-to-land ) [ 
    set disttocoast-choice max disttocoast-ahead 
    set i 0
    while [ i < length disttocoast-ahead ] [ 
      if ( item i disttocoast-ahead = disttocoast-choice ) [ set sel-angle item i angles ]
      set i i + 1
    ]
  ]
  set heading heading + sel-angle
;  if (debug = 7 and (who = 0 or who = 1) ) [ 
;    type word "Heading C: " heading 
;    print word ", angle: " sel-angle 
;  ]


  
  ; ### when to stop dispersing:  ###
  
  ; shift to dispersal away from here, along coast, if porp cannot move for one day  
  if ( distancexy (item 0 (item 1 pos-list-daily)) (item 1 (item 1 pos-list-daily)) < 2 ) [ 
    set disp-type 2
    if (debug = 7 and who = 0 ) [ 
      print "Not mov 1 d, chg to disp-type 2 (porp 0)"  
    ]
    if (debug = 7 and who = 1 ) [ 
      print "Not mov 1 d, chg to disp-type 2 (porp 1)"  
    ]
  ]

  ; shift to dispersal away from here, along coast, if porp moves too little
  if ( distancexy (item 0 (item 8 pos-list-daily)) (item 1 (item 8 pos-list-daily)) < 6 ) [ 
    set disp-type 2
    if (debug = 7 and (who = 0 or  who = 1) ) [ 
      show "Mov <15 km in 8d, chg to disp-type 2 (porp 1)"  
    ]
  ]

  ; Close to coast, chg. disp. mode 
  if ( [ disttocoast ] of patch-here < min-dist-to-land ) [ 
    set disp-type 2
    if (debug = 7 and (who = 0 or who = 1 )) [ 
      show "TOO CLOSE to land, chg to disp-type 2"  
    ]
  ]

  ; Close to target, chg. disp. mode 
  if ( the-dist < 50 ) [   ; each block is 100 x 100 cells
    set disp-type 2
    if (debug = 7 and (who = 0 or who = 1) ) [ 
      show "Close to target, chg to disp-type 2"  
    ]
  ]     
  
  if ( not ( [ bathymetry ] of patch-ahead (mean-disp-dist * 4) > 0 and [ bathymetry ] of patch-ahead (mean-disp-dist * 3) > 0 and [ bathymetry ] of patch-ahead (mean-disp-dist * 2) > 0 ) ) [    
    set disp-type 2
    if (debug = 7 and (who = 0 or who = 1) ) [ 
      show "LAND ahead, chg to disp-type 2 (porp 1)"  
    ]
  ]

  if ( not enough-water-ahead ) [
    set disp-type 0
    if (debug = 7 and who = 0 ) [ 
      print "NO WATER, stop dispersing (1) (porp 0)"
      pen-up
    ]
    if (debug = 7 and who = 1 ) [ 
      print "NO WATER, stop dispersing (1) (porp 1)"
      pen-up
    ]
  ]
  
  if ( disp-type = 1 ) [ 
    fd ( mean-disp-dist / 0.4 )
    if not ([bathymetry] of patch-here >= 0) [ ; catching rare errors
      fd (- mean-disp-dist / 0.4)
      set disp-type 0
    ]
    set energy-level energy-level - 0.001 * e-use-per-km * mean-disp-dist / 0.4   ; distances are measured in number of cells here
  ]  

end ; end porp-disp-1


to porp-disp-2  ; disperse along coast, away from prev position
  ; Disperse at least this dist from coast (same var name as in porp-disp-1)
  let min-dist-to-land 500 

  if (debug = 7 and (who = 0) ) [ 
    set color blue
    ask ( porp 0 ) [ pen-down ] 
  ]
  if (debug = 7 and (who = 1) ) [ 
    set color blue
    ask ( porp 1 ) [ pen-down ] 
  ]
  
  ; Turn away from place visited 1 day ago
  facexy (item 0 (item 1 pos-list-daily)) (item 1 (item 1 pos-list-daily))
  rt 180

  ; adjust angle to swim at const dist from land
  let disttocoast-ahead ( list
    [ disttocoast ] of patch-at-heading-and-distance (heading - 80) mean-disp-dist
    [ disttocoast ] of patch-at-heading-and-distance (heading - 70) mean-disp-dist
    [ disttocoast ] of patch-at-heading-and-distance (heading - 60) mean-disp-dist
    [ disttocoast ] of patch-at-heading-and-distance (heading - 50) mean-disp-dist
    [ disttocoast ] of patch-at-heading-and-distance (heading - 40) mean-disp-dist
    [ disttocoast ] of patch-at-heading-and-distance (heading - 30) mean-disp-dist
    [ disttocoast ] of patch-at-heading-and-distance (heading - 20) mean-disp-dist
    [ disttocoast ] of patch-at-heading-and-distance (heading - 10) mean-disp-dist
    [ disttocoast ] of patch-at-heading-and-distance (heading) mean-disp-dist
    [ disttocoast ] of patch-at-heading-and-distance (heading + 10) mean-disp-dist
    [ disttocoast ] of patch-at-heading-and-distance (heading + 20) mean-disp-dist
    [ disttocoast ] of patch-at-heading-and-distance (heading + 30) mean-disp-dist
    [ disttocoast ] of patch-at-heading-and-distance (heading + 40) mean-disp-dist
    [ disttocoast ] of patch-at-heading-and-distance (heading + 50) mean-disp-dist
    [ disttocoast ] of patch-at-heading-and-distance (heading + 60) mean-disp-dist
    [ disttocoast ] of patch-at-heading-and-distance (heading + 70) mean-disp-dist
    [ disttocoast ] of patch-at-heading-and-distance (heading + 80) mean-disp-dist
  )
  
  
;  let low-dtc-ahead not (item 1 disttocoast-ahead > 500 and item 2 disttocoast-ahead > 500 and item 3 disttocoast-ahead > 500 and item 4 disttocoast-ahead > 500 and item 5 disttocoast-ahead > 500) ; that is, if some pos ahead are close to land
;  let old-heading heading
;  facexy (item 0 (item 2 pos-list-daily)) (item 1 (item 2 pos-list-daily))
;  let tmp-heading heading
;  set heading old-heading
;  let tmp-angle subtract-headings tmp-heading old-heading  ; how much should I turn to go towards place I visited 2 days ago?
;  ifelse (tmp-angle < 0 ) [ rt -90 ] [ rt 90 ] ; turn left if that gets you back where you were previously
 
  
;  if (debug = 7 and (who = 0 or who = 1) ) [ 
;    print ""
;    print word "disttocoast-here: " [ disttocoast ] of patch-here 
;    print word "disttocoast-ahead: " disttocoast-ahead 
;  ]
  
  ; Stay on current dist from land if 1-4 km from land, or try to get there
  let dtc-diff 9999999 ; (list 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 )
  let dtc-max -9999
  let dtc-min 9999
  let dtc-nbr-alike 9
  let dtc-nbr-max 9
  let dtc-nbr-min 9
  let tmp 99999
  let dtc-here [ disttocoast ] of patch-here

  let i 0
  while [ i < length ( disttocoast-ahead ) ] [
    set tmp abs( item i disttocoast-ahead  - dtc-here )
    if (tmp <= dtc-diff) [ ; if turning this much gets you to stay on same depth, choose that angle number
      set dtc-diff tmp
      set dtc-nbr-alike i
    ]
    ; set tmp item i disttocoast-ahead
    if (tmp <= dtc-min) [
      set dtc-min tmp 
      set dtc-nbr-min i
    ]
    if (tmp > dtc-max) [ 
      set dtc-max tmp 
      set dtc-nbr-max i
    ]

    set i i + 1
  ]
  
  let angles (list -80 -70 -60 -50 -40 -30 -20 -10 0 10 20 30 40 50 60 70 80)  
  let sel-angle -9999
  if ( [disttocoast] of patch-here > 4000 ) [ set sel-angle item dtc-nbr-min angles ]
  if ( [disttocoast] of patch-here < 1000 ) [ set sel-angle item dtc-nbr-max angles ]
  if ( [disttocoast] of patch-here <= 4000 and [disttocoast] of patch-here >= 2000) [ set sel-angle item dtc-nbr-alike angles ]
  set heading heading + sel-angle

;  if (debug = 7 and (who = 0 or who = 1) ) [ 
;    print word "dtc-here: " dtc-here
;    ; print disttocoast-ahead
;    ; print word "angles: " angles
;    print word "selected angle: " sel-angle 
;    ; print ""
;  ]

  
  ; ### when to stop dispersing:  ###
    
  if ( not enough-water-ahead ) [    ; distances are measured in number of cells here
    set disp-type 0
    if (debug = 7 and who = 0 ) [ 
      print "NO WATER, stop dispersing (2 - not enough water ahead) (porp 0)"
      pen-up
      set color red
    ]
    if (debug = 7 and who = 1 ) [ 
      print "NO WATER, stop dispersing (2 - not enough water ahead) (porp 1)"
      pen-up
      set color orange
    ]
  ]


  if ( not ( [ bathymetry ] of patch-ahead ( mean-disp-dist / 0.4 ) > min-depth ) ) [    ; distances are measured in number of cells here
    set disp-type 0
    if (debug = 7 and who = 0 ) [ 
      show "LOW water, stop dispersing (2) "
      pen-up
      set color red
    ]
    if (debug = 7 and who = 1 ) [ 
      show "LOW water, stop dispersing (2)"
      pen-up
      set color orange
    ]
  ]


  if ( disp-type = 2 ) [ 
    fd ( mean-disp-dist / 0.4 )
    set energy-level energy-level - 0.001 * e-use-per-km * mean-disp-dist / 0.4
  ]     
end ; end porp-disp-2


to porp-upd-pregnancy-status
  ; 0 (unable to mate, young/low energy); 1 (unable to mate, pregnant); 2 (ready to mate)
  ;let birth-rate (list 0 0 0.136 0.417 0.818 0.714 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833 0.833) 
  ; Birth-rate data from Read, A. J. 1990. Age at sexual maturity and pregnancy rates of harbour porpoises, Phocoena phocoena, from the Bay of Fundy. Canadian Journal of Fisheries and Aquatic Sciences 47:561-565.

  ; Become ready to mate:
  if (pregnancy-status = 0 and age >= age-of-maturity ) [ set pregnancy-status 2 ]

  ; Mate:
  if (pregnancy-status = 2 and round (sim-day - 360 * (year - 1)) = mating-day) [ 
    if (random-float 1 < 0.68) [ ;  become pregnanat with prob. taken from Read & Hohn 1995
      set pregnancy-status 1
      if (debug = 9 ) [ print word who " pregnant" ]
      set ds-mating 0
    ]
  ]
  
  ; Give birth:
  if (pregnancy-status = 1 and ds-mating = 10 * 30) [ ; give birth. Gestation time = approx 10 mo (Lockyer 2003)
    set pregnancy-status 2  ; so it is ready to mate even though it has a very young calf
    set with-lact-calf true
    set ds-mating -99
    set dsg-birth 0
    if (debug = 9 ) [ print word who " with lact calf" ]
  ]
  
  if (with-lact-calf = true and dsg-birth = 8 * 30) [  ; nursing for 8 months
    set with-lact-calf false
    set dsg-birth -99
    let n-offspr 0
    if (random 2 = 1 ) [ set n-offspr 1 ]  ; assuming 50 % males and no abortions
    if (debug = 9 ) [ 
      let tmp word who " hatching "
      print word tmp n-offspr
    ]
    hatch-porps  n-offspr [  
        set age 0
        set pregnancy-status 0
        set energy-level random-normal 10 1
    ]
  ]

  if (pregnancy-status = 1) [ set ds-mating ds-mating + 1 ]
  if (with-lact-calf) [ set dsg-birth dsg-birth + 1 ]
end ; end porp-upd-pregnancy-status

to porp-upd-mortality
  ; Introducing maximum age
  if (age > 30) [ 
    set list-of-dead-age lput (floor age) list-of-dead-age
    set list-of-dead-day lput (floor sim-day) list-of-dead-day
    die 
  ]
  
  ; Mortality due to by-catch
  let daily-survival-prob exp( ln((1 - bycatch-prob)) / 360 )  ; Ok that only divided by 360, called once per day
  if (random-float 1 > daily-survival-prob) [
    set list-of-dead-age lput (floor age) list-of-dead-age
    set list-of-dead-day lput (floor sim-day) list-of-dead-day
    die    
  ]
  
end ; end porp-upd-mortality


; Statistics and plots
 
to my-update-plots ; update histograms and other plots 
  ; food level and number of porpoises
  set-current-plot "population"
  set-current-plot-pen "N x10"
  plot 10 * ( count porps )
  set-current-plot-pen "total food"
  plot sum [ food-level ] of patches with [ food-level > 0 ]
  set-current-plot-pen "E x100"
  plot 100 * mean [ energy-level ] of porps

  ; Histogram of energy-level in porpoises
  set-current-plot "porpoise-energy"
  histogram [ energy-level ] of porps
  
  ; Age-class distribution in porpoise population
  set-current-plot "age-distribution"
  histogram [ age ] of porps

end


; File I/O

to file-setup
  let repl-counter 0                            ; different replicates of the same porpoise
  let go-on true
  while [go-on] [
    set repl-counter repl-counter + 1
    let repl-counter-2 word "1000" repl-counter
    set repl-counter-2 substring repl-counter-2 (length repl-counter-2 - 3) (length repl-counter-2)
    let file-tmp word "output/" area
    set file-tmp word file-tmp "/"
    set file-tmp word file-tmp output-name
    set file-tmp word file-tmp "-"
    set file-tmp word file-tmp repl-counter-2
    set outfile word file-tmp ".txt"
    if (not file-exists? outfile) [ 
      set go-on false
      file-open outfile
    ]
    if (repl-counter = 30) [ 
      set go-on false 
      user-message ( "The desired number of replicates has been obtained" )
    ]
  ]
  file-print (" id day year utm-x utm-y block age energy") ; header line (space-separated). Note that "bathy" is non-standard.
end

to file-setup-2
  let repl-counter 0                            ; different replicates of the same porpoise
  let go-on true
  while [go-on] [
    set repl-counter repl-counter + 1
    let repl-counter-2 word "1000" repl-counter
    set repl-counter-2 substring repl-counter-2 (length repl-counter-2 - 3) (length repl-counter-2)
    let file-tmp word "output/" area
    set file-tmp word file-tmp "/"
    set file-tmp word file-tmp output-name
    set file-tmp word file-tmp "-"
    set file-tmp word file-tmp repl-counter-2
    set outfile word file-tmp ".txt"
    if (not file-exists? outfile) [ 
      set go-on false
      file-open outfile
    ]
    if (repl-counter = 30) [ 
      set go-on false 
      user-message ( "The desired number of replicates has been obtained" )
    ]
  ]
  file-print (" day utm-x utm-y block") ; header line (space-separated). Note that "bathy" is non-standard
end


to file-write-line
  ; Ask porpoise to write variables that should be imported into an R track object:
  ; "animal ptt pop sex length weight x y year month day hour minute second":
  file-open outfile ; append
  ; add entries
  ask porps [
    file-write who ; not followed by CR
    file-write round sim-day ; not followed by CR
    file-write year
    file-write round(xcor * 400 + xllcorner) ; 240 km wide and 400 km tall non-wrapped landscape divided into 400 x 400 m cells
    file-write round(ycor * 400 + yllcorner)
    file-write [ block ] of patch-here 
    ;file-write lgth
    ;file-write weight
    file-write age    
    file-write energy-level
    file-print ""    ; followed by CR
  ]
  file-close
end


to file-write-line-2
  ; Ask porpoise to write variables that should be imported into an R track object:
  ; "animal ptt pop sex length weight x y year month day hour minute second":
  file-open outfile ; append
  ; add entries
  ask porp 0 [
    file-write round sim-day ; not followed by CR
    file-write round(xcor * 400 + xllcorner) ; 240 km wide and 400 km tall non-wrapped landscape divided into 400 x 400 m cells
    file-write round(ycor * 400 + yllcorner)
    file-write [ block ] of patch-here 
;    file-write lgth
;    file-write weight
;    file-write age    
;    file-write energy-level
    file-print ""    ; followed by CR
  ]
  file-close
end


; Main

to setup
  clear-turtles
  clear-output
  clear-drawing   ; clears pendown tracks
  clear-all-plots
  reset-ticks
  landsc-setup
  porps-setup
  update-block-values
  if ( model >= 3 ) [ 
    if (not (wind-farms = "off")) [
      turbs-import-pos 
      turbs-setup
    ]
    if (incl-ships) [
      ships-import
      ships-setup 
    ]
  ]
  if ( write-data != "off" and write-data != "one-porp") [ 
    file-setup
    print word "Writing to: " outfile
  ]
  
;  ; Choose daily survival probabilities
;  ; Max, mean and min curves for yearly survival for 9 representative species from Caswell 1998 converted to daily survival values (not energy related)
;  ; Each value in the lists is the daily mean used for porps in a particular age class (same val used for a whole year)
;  if (survival = "high") [ set daily-surv-prob (list 0.999770956 0.999814739 0.999835458 0.999931317 0.999965003 1 1 1 0.999964552 0.999964088 0.999963613 0.999925748 0.99988472 0.99987967 0.999830899 0.999772851 0.999752363 0.999727809 0.999697845 0.999660459 0.999527995 0.999322043 0.999387156 0.999387156 0.999429725 0.999279479 0.999020604 0.999076053 0.997485883 0.998097561 0) ]
;  if (survival = "mean") [ set daily-surv-prob (list 0.998358939 0.999509209 0.999546587 0.999541166 0.99954872 0.999576599 0.999633223 0.999576599 0.999676473 0.999633223 0.999576599 0.999499242 0.99971059 0.999309814 0.999576599 0.999499242 0.999387156 0.999209977 0.998886705 0.999209977 0.998886705 1 0.998097561 0 0 0 0 0 0 0 0) ]
;  if (survival = "low") [ set daily-surv-prob (list 0.99915452 0.999704355 0.999738193 0.999810477 0.999848748 0.999894329 0.99991799 0.999915466 0.999912782 0.999909922 0.999906868 0.999870702 0.999864315 0.999820393 0.999807825 0.999750109 0.999725083 0.999694482 0.999656204 0.999641927 0.999506185 0.999345239 0.999347269 0.999319363 0.999429725 0.999279479 0.999020604 0.999076053 0.997485883 0.998097561 0) ]

  
  ; Debugging code:
  if ( write-data = "one-porp") [ 
    file-setup-2
    print word "Writing porp 0 data to: " outfile
  ]
  if (debug = 1) [ 
    print ""
    print "Debug avoidance behaviour near land (debug=1):"
  ]
  if (debug = 2) [ 
    print ""
    print "Write turning angles before/after react. to land (debug=2):"
  ]
  if (debug = 3) [ 
    print ""
    print "Debugging CRW related turns (mod >=1) (debug=3):"
  ]
  if ( (debug = 4) and (model >= 2) ) [ 
    print ""
    print "Debugging attraction vector (mod >=2) (debug=4):"
  ]
  if ( (debug = 5) and (model >= 2) ) [ 
    print ""
    print "Debugging attraction + deterrence vectors (mod >=2)  (debug=6):"
  ]
  if ( debug = 6 ) [ 
    print "" 
    print "Showing mean energy per 40 x 40 km block (debug=6)"
    print "day \tSE Anh (16) \tGB-N (33)  \tNV Born (42) \tFehm (52)"
  ]
  if (debug = 7) [ 
    print "Debugging dispersal behaviour" 
    inspect porp 0    
    inspect porp 1  
  ]
  if (debug = 8) [ print "Debugging deterrence from wind farms and ships" ]
  if (debug = 9) [ 
    print "Debugging reprod. and mortality" 
  ]
end


to daily-tasks
    ; Things to do daily or less often
    my-update-plots
    landsc-display
    
    ; print word "age of dead: " list-of-dead-age
    
    ; write data to file:
    if write-data = "daily" [ ; must write deployment pos before moving in movement models
      ask ( porps ) [ file-write-line ] 
    ]
    if write-data = "one-porp" [ ; one porpoise (porp 0), daily
      ask ( porp 0 ) [ file-write-line-2 ] 
    ]
    if ( not (prev-year = year ) or ( time-step = 1 ) ) [
      if write-data = "yearly" [
        ask ( porps ) [ file-write-line ] 
      ]
    ]

    
    ; Update daily average energy level and corresponding positions for porps and use it to start/stop dispersing
    ask ( porps ) [
      set age age + (1 / 360)
      let e-list-lgt length energy-level-daily 
      set energy-level-daily remove-item ( e-list-lgt - 1 ) energy-level-daily
      let e-mean energy-level-sum / 48  ; 48 half-hour steps per day
      set energy-level-daily fput (precision e-mean 3) energy-level-daily
      set pos-list-daily fput ( list xcor ycor ) pos-list-daily
      set pos-list-daily remove-item ( length pos-list-daily - 1) pos-list-daily
      
      if ( model >= 3 ) [
        if ( disp-type = 0 ) [ ; i.e. not dispersing, 
          if ( ( item 0 energy-level-daily ) < ( item 1 energy-level-daily ) and ( item 1 energy-level-daily ) < ( item 2 energy-level-daily )  and ( item 2 energy-level-daily ) < ( item 3 energy-level-daily ) ) [ 
            set disp-type 1 ; decreasing energy for three days
            porp-disp-target-select
          ]
        ]
        if ( disp-type > 0 ) [
          ; Energy level higher than any of the previous seven days, stop dispersing;
          if ( item 0 energy-level-daily ) > max (list item 1 energy-level-daily item 2 energy-level-daily item 3 energy-level-daily item 4 energy-level-daily item 5 energy-level-daily item 6 energy-level-daily item 7 energy-level-daily ) [ 
            set disp-type 0 
            if (debug = 7 and who = 0) [ 
              ask ( porp 0 ) [ pen-up ]
              set color red
              print  "Food found, stop disp., (porp 0)"
            ]
            if (debug = 7 and who = 1) [ 
              ask ( porp 1 ) [ pen-up ]
              set color red
              print  "Food found, stop disp., (porp 1)"
            ]
            set disp-target list 0 0
          ]
        ]
        if (disp-type = 2) [
          ; Energy level higher last week, turn towards place visited 3 days ago
          if ( mean( list(item 0 energy-level-daily) (item 1 energy-level-daily) (item 3 energy-level-daily)) < mean( list(item 6 energy-level-daily) (item 7 energy-level-daily) (item 8 energy-level-daily) (item 9 energy-level-daily) ) ) [
            facexy (item 0 (item 2 pos-list-daily)) (item 1 (item 2 pos-list-daily))
            ;set heading ( heading + random 16 - 8 )
            ;rt 180
            fd 1
            if (debug = 7 and (who = 0 or who = 1)) [ show "RETURNING to prev pos" ]           
          ]
        ]
      ] ; end model >=3 tasks
      
    ]
    ask ( porps ) [
      set energy-level-sum 0 ; reset daily
    ]
    
    ; Update parameters related to demography
    if ( model >= 4 ) [
      ask ( porps ) [
        porp-upd-mortality
        porp-upd-pregnancy-status
      ]
    ]
        
    if ( debug = 6 ) [ 
      ; get mean energy of porps in blocks 16 (SE of Anholt), 33 (N Great Belt), 42 (NCV of Bornholm) and 52 (Fehmarn)
      type word round(sim-day) "\t"
      let e-list [ energy-level ] of porps with [ block = 16 ]
      if ( length e-list > 0 ) [ type word (precision mean e-list 2) "      \t" ]
      if ( length e-list = 0 ) [ type word "- " "\t\t" ]
      set e-list [ energy-level ] of porps with [ block = 33 ]
      if (length e-list > 0) [ type word (precision mean e-list 2) "      \t" ]
      if (length e-list = 0) [ type word "- " "\t\t" ]
      set e-list [ energy-level ] of porps with [ block = 42 ]
      if (length e-list > 0) [ type word (precision mean e-list 2) "      \t" ]
      if (length e-list = 0) [ type word "- " "\t\t" ]
      set e-list [ energy-level ] of porps with [ block = 52 ]
      if (length e-list > 0) [ print word (precision mean e-list 2) "      \t" ]
      if (length e-list = 0) [ print word "- " "\t\t" ]
    ]
end ; end daily-tasks

to calc-age-distrib
  ; age list at start of year
  set age-distrib [ ]
  let aa [round age] of porps
  set age-distrib lput ( length filter [? = 0] aa ) age-distrib
  set age-distrib lput ( length filter [? = 1] aa ) age-distrib
  set age-distrib lput ( length filter [? = 2] aa ) age-distrib
  set age-distrib lput ( length filter [? = 3] aa ) age-distrib
  set age-distrib lput ( length filter [? = 4] aa ) age-distrib
  set age-distrib lput ( length filter [? = 5] aa ) age-distrib
  set age-distrib lput ( length filter [? = 6] aa ) age-distrib
  set age-distrib lput ( length filter [? = 7] aa ) age-distrib
  set age-distrib lput ( length filter [? = 8] aa ) age-distrib
  set age-distrib lput ( length filter [? = 9] aa ) age-distrib
  set age-distrib lput ( length filter [? = 10] aa ) age-distrib
  set age-distrib lput ( length filter [? = 11] aa ) age-distrib
  set age-distrib lput ( length filter [? = 12] aa ) age-distrib
  set age-distrib lput ( length filter [? = 13] aa ) age-distrib
  set age-distrib lput ( length filter [? = 14] aa ) age-distrib
  set age-distrib lput ( length filter [? = 15] aa ) age-distrib
  set age-distrib lput ( length filter [? = 16] aa ) age-distrib
  set age-distrib lput ( length filter [? = 17] aa ) age-distrib
  set age-distrib lput ( length filter [? = 18] aa ) age-distrib
  set age-distrib lput ( length filter [? = 19] aa ) age-distrib
  set age-distrib lput ( length filter [? = 20] aa ) age-distrib
  set age-distrib lput ( length filter [? = 21] aa ) age-distrib
  set age-distrib lput ( length filter [? = 22] aa ) age-distrib
  set age-distrib lput ( length filter [? = 23] aa ) age-distrib
  set age-distrib lput ( length filter [? = 24] aa ) age-distrib
  set age-distrib lput ( length filter [? >= 25] aa ) age-distrib
end ; end calc-age-distrib

to calc-mort-prob
  ; age distribution in list of dead
  let mort-age-distrib [ ]
  set mort-age-distrib lput ( length filter [? = 0] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [? = 1] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [? = 2] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [? = 3] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [? = 4] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [? = 5] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [? = 6] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [? = 7] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [? = 8] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [? = 9] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [? = 10] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [? = 11] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [? = 12] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [? = 13] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [? = 14] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [? = 15] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [? = 16] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [? = 17] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [? = 18] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [? = 19] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [? = 20] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [? = 21] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [? = 22] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [? = 23] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [? = 24] list-of-dead-age ) mort-age-distrib
  set mort-age-distrib lput ( length filter [? >= 25] list-of-dead-age ) mort-age-distrib
  
  ; report results
  if (debug = 0 and model >= 4) [
    if (sim-day < 2) [ 
      print " "
      print "Age class distrib. and yearly mort.:"
      print " y age n mort"
    ]
    foreach [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25] [
      write round year - 1
      write ?
      write (item ? age-distrib)
      print word " " (item ? mort-age-distrib) 
    ]
  ]
  
  ; reset mort lists every year, after reporting:
  set list-of-dead-age [ ]
  set list-of-dead-day [ ]
end ; end calc-mort-prob

to yearly-tasks
  ask ( porps ) [ set mating-day round( random-normal (7.5 * 360 / 12) 20 )] 
  ; Make list of number of dead per age class and corresp pop nos in prev year
  if (sim-day < 2) [ calc-age-distrib ]
  calc-mort-prob
  calc-age-distrib
end ; end yearly-tasks


to go
  if ( debug = "profile" ) [ 
    profiler:start
  ]
  if ( time-step = 0 ) [ 
    reset-timer
  ]

  if (model >= 3) [  ; deterrence behaviour
    porps-upd-deter
    turbs-deter-porps
    ships-deter-porps
  ]
  
  ask porps [
    porp-move  ; This is the important step!
    if ( disp-type = 1 ) [ porp-disp-1 ] ; long-dist dispersal, set in daily-tasks
    if ( disp-type = 2 ) [ porp-disp-2 ] ; dispersal away from old pos, along coast
    if (debug = 8) [     ; if no longer scared by turbines, change color to std
      if (deter-strength = 0) [ ifelse(who = 0) [set color red] [set color orange]] 
    ]
  ]
  
  if ( model >= 3 and incl-ships) [ ships-move ]

  if (debug = 8) [ ; write porp xy near wind farms
    if (area = "Homogeneous" and wind-farms = "Line" and write-data = "off") [ file-noise-debug ]
  ]
  
  set time-step time-step + 1
  let prev-day floor sim-day
  set prev-month month
  set prev-quarter quarter
  set prev-year year
  set sim-day time-step / 48                                          ; day since start of sim
  set year floor ( sim-day / 360 ) + 1                                ; this gives 90 days per quarter
  set month ceiling ( (sim-day - (year - 1) * 360) / 30 ) 
  set quarter 1 + floor ( (sim-day + 30 - (year - 1) * 360) / 90 )    ; quarter 1 is Dec-Febr etc.
  if quarter = 5 [ set quarter 1 ]
  ;tick                                                               ; slows things down a lot

  ; Update food growth map and food level
  if ( not (prev-quarter = quarter )  ) [
    landsc-upd-maxent-level-map
  ]

  if ( remainder sim-day food-upd-interval ) = 0 [
    landsc-upd-food
  ]

  ; Things to do yearly: 
  if ( not (prev-year = floor year ) or ( time-step = 1 ) ) [
    yearly-tasks
  ]

  ; Things to do daily (or less often): 
  if ( not (prev-day = floor sim-day ) or ( time-step = 1 ) ) [
    daily-tasks
  ]

  if ( not (prev-month = month ) or ( time-step = 1 ) ) [
    if (write-data = "monthly") [
      ask ( porps ) [ file-write-line ]
    ]
  ]


  ; Stop running and print summary
  if ( time-step / 48 = max-sim-day ) [
    let ttt word "Time (" time-step
    set ttt word ttt " half-hour intervals): "
    set ttt word ttt timer
    print word ttt " sec"
    ; print ( list "Time (" max-tick "ticks): " timer "sec" )
    if (write-data != "off") [ file-close  ]
    stop
  ]

  if ( debug = "profile" ) [ 
    profiler:stop
    print profiler:report
    profiler:reset
  ]
end
@#$#@#$#@
GRAPHICS-WINDOW
382
11
992
1042
-1
-1
0.5
1
7
1
1
1
0
1
1
1
0
599
0
999
0
0
0
ticks
30.0

BUTTON
13
207
148
240
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
14
248
77
281
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

CHOOSER
158
196
298
241
disp-var
disp-var
"bathymetry" "maxent-level" "food-prob" "food-level" "blocks"
3

MONITOR
305
11
368
56
NIL
time-step
0
1
11

BUTTON
85
248
148
281
NIL
go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

MONITOR
305
67
368
112
NIL
sim-day
1
1
11

CHOOSER
118
312
210
357
debug
debug
"profile" 0 1 2 3 4 5 6 7 8 9
1

CHOOSER
16
312
108
357
model
model
0 1 2 3 4
4

SLIDER
13
11
148
44
n-porps
n-porps
1
500
200
1
1
NIL
HORIZONTAL

CHOOSER
13
52
148
97
area
area
"Kattegat" "Homogeneous"
0

BUTTON
158
248
298
281
Update disp-var
landsc-display
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
21
934
191
967
food-growth-rate
food-growth-rate
0
0.2
0.1
0.01
1
(rU)
HORIZONTAL

MONITOR
305
124
368
169
NIL
year
17
1
11

MONITOR
305
180
368
225
NIL
quarter
1
1
11

INPUTBOX
158
11
298
71
max-sim-day
14400
1
0
Number

PLOT
17
369
371
553
population
time (days)
NIL
0.0
900.0
0.0
4000.0
true
true
"" ""
PENS
"N x10" 1.0 0 -10899396 true "" ""
"E x100" 1.0 0 -7500403 true "" ""
"total food" 1.0 0 -2674135 true "" ""

MONITOR
306
501
365
546
N porps
count porps
17
1
11

PLOT
18
566
196
716
porpoise-energy
NIL
NIL
0.0
20.0
0.0
80.0
true
false
"" ""
PENS
"default" 1.0 1 -16777216 true "" ""

TEXTBOX
22
915
186
933
Food growth and energy use
11
0.0
1

SLIDER
21
973
191
1006
e-use-per-km
e-use-per-km
0
1
0
0.1
1
E-3
HORIZONTAL

SLIDER
21
1012
191
1045
e-use-per-30min
e-use-per-30min
4
8
4.5
0.1
1
E-3
HORIZONTAL

CHOOSER
158
143
298
188
write-data
write-data
"off" "daily" "monthly" "yearly" "one-porp"
0

INPUTBOX
158
77
298
137
output-name
kat-111220
1
0
String

SLIDER
19
829
190
862
mean-disp-dist
mean-disp-dist
0.1
6
1.6
0.1
1
km/30min
HORIZONTAL

TEXTBOX
25
732
175
750
Dispersal parameters
11
0.0
1

SLIDER
19
868
190
901
min-dist-to-target
min-dist-to-target
20
200
100
20
1
km
HORIZONTAL

SLIDER
19
750
190
783
min-disp-depth
min-disp-depth
0
10
4
0.1
1
m
HORIZONTAL

TEXTBOX
19
293
169
311
Model selection and testing
11
0.0
1

CHOOSER
13
104
149
149
wind-farms
wind-farms
"off" "Pot_Krieger" "Pot_St_Middelgr" "Rodsand-I" "Rodsand-II" "Samsoe" "Sprogoe" "User-def" "All" "Line" "Potential"
0

SWITCH
13
156
148
189
incl-ships
incl-ships
0
1
-1000

SLIDER
200
895
371
928
deterrence-coeff
deterrence-coeff
0
10
7.6
0.1
1
NIL
HORIZONTAL

TEXTBOX
201
877
351
895
Noise avoidance param.
11
0.0
1

SLIDER
200
934
371
967
std-deterrence-dist
std-deterrence-dist
200
1600
300
20
1
m
HORIZONTAL

SLIDER
200
974
371
1007
deter-time
deter-time
1
10
5
1
1
time steps
HORIZONTAL

MONITOR
306
452
365
497
N lact
count porps with [ with-lact-calf = true ]
0
1
11

PLOT
202
566
371
716
age-distribution
year class
NIL
0.0
25.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 1 -16777216 true "" ""

MONITOR
305
236
368
281
NIL
month
1
1
11

TEXTBOX
202
732
352
750
Reproduction and mortality
11
0.0
1

TEXTBOX
312
433
371
451
in patches
11
0.0
1

SLIDER
19
789
190
822
n-disp-targets
n-disp-targets
4
12
12
2
1
NIL
HORIZONTAL

SLIDER
200
750
371
783
x-survival-prob-const
x-survival-prob-const
0.1
2
0.4
0.1
1
NIL
HORIZONTAL

SLIDER
200
829
371
862
bycatch-prob
bycatch-prob
0
0.2
0
0.001
1
y-1
HORIZONTAL

SLIDER
200
789
371
822
m-mort-prob-const
m-mort-prob-const
0
1
1
0.01
1
NIL
HORIZONTAL

SWITCH
200
1013
371
1046
pile-driving
pile-driving
1
1
-1000

@#$#@#$#@
## ODD FOR THE HARBOUR PORPOISE POPULATION MODEL

Model and description updated by Jacob Nabe-Nielsen 2011-04-05

## PURPOSE

The purpose of the model is to investigate the relative impact of wind farms and ship traffic on the dynamics of Danish harbour porpoise populations. The noise from wind farms and ships may scare porpoises away and thereby cause habitat fragmentation and reduced amounts of available food, which is likely to affect porpoise population sizes.

## STATE VARIABLES AND SCALES

Individuals are characterised by the state variables: age, energetic status, pregnancy status and lactating. Each individual in the model represent 20 female porpoises. Animals in the age class 8 months to 3.44 years are independent juveniles (cf. Read 1990 and Lockyer and Kinze. 2003). Younger animals (calves) are not included as independent individuals. The pregnancy status of independent animals can be pregnant/not pregnant, or infertile. Individuals are assumed not to interact except through their consumption of a common resource.

Anthropogenic objects (AOs) are characterised by the state variables: type and noise level. Two types of AO are modelled: wind turbines and boats (that are able to move).

Simulations are based on a 240 km wide and 400 km tall non-wrapped landscape divided into 400 x 400 m cells (i.e. 600x1000 cells cells in total) and sixty 40 x 40 km blocks. The landscape represents Kattegat and the waters between Denmark and Germany. It includes three kinds of environment: land (52.1%), water without food, and 1-cell large food patches. The amount of food in a patch is characterised by a variable food level.

## PROCESS OVERVIEW AND SCHEDULING

The model runs in half-hour steps, and individuals respond to land and AOs by turning after each step. Fine-scale animal movement is modelled to result from a mixture of correlated random walk (CRW) behaviour and a memory-based ability to return to areas where individuals found food previously (see Nabe-Nielsen et al., submitted). 

The animal energy level E is updated after each time step. The individuals� energy consumption per half-hour step is modelled to reflect basic metabolic costs, which are higher when the water is cold (approx. 15% increase in April and October and 30% in May�September; Lockyer et al. [2003]) and to increase when animals lactate (40% increase, M. Wahlberg, pers. comm.). The individuals� energy E level is scaled to lie in the range 0�20 in the model. When their energy level is higher than 10 they build up energy reserves, causing them to consume a smaller proportion of the food they encounter. This proportion is modelled to decrease linearly as the energy level increases from 10 to 20. The amount of food in the food patches is adjusted accordingly. Afterwards the food level increases logistically up to a maximum level; see details under Input. 

Animal mortality is modelled to depend on their energetic status, with a yearly survival probability equal to 1 � exp(�k x E), where k is a positive constant. Animals risk dying each time step after updating their energetic status, and also die when reaching 30 years. Lactating animals abandon their offspring rather instead of dying, unless their energy level is �0. Further, the model includes an age and energy independent by-catch rate of 1.7% (the Ascobans safety limit for by-catch).

Animal reproduction is divided in three phases: mating, gestation and nursing. Mating peaks in August (Lockyer 2003), and in our model each individual has a mating day that is selected as a random normal variable with mean 7.5 x 30 and a standard deviation of 20. Individuals that are sexually mature, i.e. at least 3.44 years old (Read 1990) they are modelled to become pregnant with a probability of 0.68, following Read & Hohn (1995). After 10 months they give birth to a calf (Lockyer 2003). After a lactation period of eight months the calf gives rise to a new, independent individual with a probability of 0.5 (assuming equal sex ratio). Abortions of unborn calves is not taken into account.

Long-distance dispersal results from two different processes in the model: (1) if the average daily energy level decreases for three consecutive days, porpoises turn towards a 40 x 40 km block selected at random among the 12 blocks with highest expected quality (based on average Maxent value for the blocks each quarter; see Edr�n et al. 2010) >40 km away. Afterwards they turn �20� in the direction with deepest water, provided that there is no land further away (8x disp-dist) in that direction. Finally they turn �30� to get as far away from land as possible if water depth <min-disp-dept or if <2 km from the coast. (2) When approaching the target block or if they are unable to get to an area with deep water (>min-disp-dept) the porpoises start moving away from the areas they visited the previous day. They attempt to stay at a constant distance from land if 1-4 km from land, else they try to get there. If the average energy level was higher 6�9 days ago than for the last three days the porpoise turns towards the place where it were three days ago. Finally the dispersing porpoises move disp-dist forward. The porpoises stop dispersing when they get trapped in areas with low water or when the current energy level is higher than at any time during the previous week. The target blocks are not selected entirely at random for animals that start dispersing immediately north of Djursland and Funen or E of Sealand, as these do not use directed dispersal.

Deterrence behaviour, i.e. the porpoises' reaction to wind turbines and ships, is related to the distance to the disturbing object (DO), to how much the object disturbs (its impact on the porpoise, IP) and to the maximum deterrence distance for a standard wind turbine, MD. Standard wind turbines (e.g. R�dsand II turbines) have IP=1. Whenever a porpoise gets closer than MD * IP from a ship or wind turbine, its movement direction in the next step is calculated as the sum of a deterrence vector and the vector defining the normal fine-scale movement. The deterrence vector is calculated as j * DC*(IP * MD - DO), where j is a unity vector pointing away from the disturbing object and DC is a constant deterrence coefficient that controls the balance between the standard fine-scale move and the deterrence effect throughout. The step length is not affected by the strength of the deterrence. The length of the deterrence vector therefore decreases linearly with distance to the disturbing object. This is similar to sound in water under some circumstances. At the end of each time step the length of the deterrence vector is halved, and after deter-time = 5 time steps it is set to zero (i.e. the porpoise only moves away from the disturbing object for deter-time time steps).

Variables describing the state of the food patches and of the individuals (except movement and energetic status) are updated daily in the following order: (1) increase food level in patches, (2) start dispersing if energy level drops, (3) die due to age-specific background mortality and loose unborn or lactating offspring (related to energetic status), (4) update pregnancy status: mate, give birth and weane lactating calves depending on time of year. Movement and energetic status of the individuals is updated in every time step.

See separate flowchart of animal-related processes in the model.

## DESIGN CONCEPTS

Emergence: Population dynamics emerge from the behaviour of the individuals and the balance between their energy expenditure (related to time, water temperature and lactation) and to their food acquisition rate.  
Adaptation: Individuals' responses to land and AOs is fixed, and adaptation is not modelled explicitly.  
Fitness: Shifts between two different types of movement behaviour is modelled to result from optimal foraging based on an evaluation of the amount of food acquired in the past.  
Sensing: Individuals are modelled to respond to the presence of land and AOs by adjusting their movement behaviour, but the response is assumed to be independent of their state.  
Interaction: Interactions among individuals are modelled implicitly through their competition for food.  
Stochasticity: Both movement (direction and distance moved when moving locally), mortality, mating date and dispersal behaviour (which area to disperse to) depend on stochastic processes.   
Collectives (groups of individuals): Animals do not interact in the model, except that juveniles are inextricably linked to their mother till they finish lactating.  
Observation (collecting data from the IBM): The positions of the individuals are sampled monthly and compared to an independent dataset. The age-class distribution and age specific mortality is sampled yearly and compared to independent datasets.

## INITIALISATION

The model is initialised by creating 300 porpoises and 9600 randomly distributed food patches. The initial energy level of the porpoises is modelled as a random normal variable with mean 10 and standard deviation one. Each patch has a size of one cell. Both the patch size and the number of patches (on the average 1000 per 100 x 100 km) correspond to the numbers used in the movement model. The patch locations remain constant among model runs. Patches that happen to be on land are subsequently removed so only 4572 food patches are retained. The simulation is initiated on 1 January, which affects the food replenishment rate.

## INPUT

The food level in the randomly distributed food patches increases logistically after being eaten. The maximum food level is calculated as the season-specific Maxent value for the patch divided by the mean season-specific Maxent value for the entire landscape. Maxent values fall in the range 0�1, with high values in areas with environmental conditions resembling the ones found in areas with a high porpoise density (Edr�n et al. 2010). The first season covers the months December�February. The rate of increase in the amount of food is kept constant (rU) = 1.02). See details in Nabe-Nielsen et al. (submitted).

## SUBMODELS

1. The mathematical �skeleton� of the model.

2. A full model description. This version has exactly the same structure as the �skeleton� (i.e., the same subtitles and equation numbers), but now each equation and parameter is verbally explained in full detail. The second version may be made available on the internet / in a report.

The amount of energy spent (ES) by a porpoise is modelled as linear function of time, mt , so ESt = C + D mt . Here C is a constant related to basic metabolism and increased energy expenditure related to fast movements when foraging and m is measured in 100-m steps.

## REFERENCES

Caswell, H., S. Brault, A. J. Read, and T. D. Smith. 1998. Harbor porpoise and fisheries: An uncertainty analysis of incidental mortality. Ecological Applications 8:1226-1238.

Edr�n, S. M. C., M. S. Wisz, J. Teilmann, R. Dietz, and J. S�derkvist. 2010. Modelling spatial patterns in harbour porpoise satellite telemetry data using maximum entropy. Ecography 33:698-708.

Lockyer, C. 2003. Harbour porpoises (Phocoena phocoena) in the North Atlantic: biological parameters. Harbour porpoises in the North Atlantic 5:143-175.

Lockyer, C., G. Desportes, K. Hansen, S. Labbert�, and U. Siebert. 2003. Monitoring growth and energy utilization of the harbour porpoise (Phocoena phocoena) in human care. Harbour porpoises in the North Atlantic 5:143-175.

Lockyer, C., and C. Kinze. 2003. Status, ecology and life history of harbour porpoise (Phocoena phocoena), in Danish waters. Pages 143-175 in T. Haug, G. Desportes, G. A. V�kingsson, and L. Witting, editors. Harbour porpoises in the North Atlantic. The North Atlantic Marine Mammal Commission, Troms�.

Nabe-Nielsen, J., J. Tougaard, J. Teilmann, and M. Forchhammer. How switching between movement strategies can lead to home ranges and optimal foraging. Submitted to Ecology.

Read, A. J. 1990. Age at sexual maturity and pregnancy rates of harbour porpoises, Phocoena phocoena, from the Bay of Fundy. Canadian Journal of Fisheries and Aquatic Sciences 47:561-565.

Read, A. J., and A. A. Hohn. 1995. Life in the fast lane: The life history of harbor porpoises from the Gulf of Maine. Marine Mammal Science 11:423-440.
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270

@#$#@#$#@
NetLogo 5.0.2
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="110210 Test disp-dist" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="52416"/>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="3.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-sim-110210&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;one-porp&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="3.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.5"/>
      <value value="1"/>
      <value value="2"/>
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;bathymetry&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="3"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110317  Deterrence debugging" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="288"/>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;hom-110308&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="3.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="0.5"/>
      <value value="0.75"/>
      <value value="1"/>
      <value value="1.5"/>
      <value value="2"/>
      <value value="4"/>
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Homogeneous&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;Line&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
      <value value="400"/>
      <value value="500"/>
      <value value="600"/>
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110322a  Deterrence debugging" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="288"/>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;hom-110308&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="3.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="0.5"/>
      <value value="0.75"/>
      <value value="1"/>
      <value value="1.5"/>
      <value value="2"/>
      <value value="4"/>
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Homogeneous&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;Line&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110322b  Deterrence debugging" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="288"/>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;hom-110308&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="3.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="16"/>
      <value value="32"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Homogeneous&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;Line&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="200"/>
      <value value="300"/>
      <value value="400"/>
      <value value="500"/>
      <value value="600"/>
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110404a Calibrate k and e-use" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-interface word random 10000 ".png"
export-output word random 10000 ".txt"</final>
    <timeLimit steps="864000"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="abandon-juv-at-e">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="added-d-juv-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="1.8"/>
      <value value="2.2"/>
      <value value="2.6"/>
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.02"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110404b Calibrate k and e-use" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-interface word random 10000 ".png"
export-output word random 10000 ".txt"</final>
    <timeLimit steps="864000"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="abandon-juv-at-e">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="added-d-juv-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="1.8"/>
      <value value="2.2"/>
      <value value="2.6"/>
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.02"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110404c Calibrate k and e-use" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-interface word random 10000 ".png"
export-output word random 10000 ".txt"</final>
    <timeLimit steps="864000"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="abandon-juv-at-e">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="added-d-juv-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="1.8"/>
      <value value="2.2"/>
      <value value="2.6"/>
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.02"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110404d Calibrate k and e-use" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-interface word random 10000 ".png"
export-output word random 10000 ".txt"</final>
    <timeLimit steps="864000"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="abandon-juv-at-e">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="added-d-juv-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="1.8"/>
      <value value="2.2"/>
      <value value="2.6"/>
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.02"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110405 Calibrate k and e-use" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-interface word random 10000 ".png"
export-output word random 10000 ".txt"</final>
    <timeLimit steps="345600"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="added-d-juv-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110406a Calibrate k and e-use" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="m-survival-prob-const">
      <value value="0.75"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.3"/>
      <value value="0.4"/>
      <value value="0.5"/>
      <value value="0.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="2.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="added-d-juv-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110406b Calibrate k and e-use" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="m-survival-prob-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.3"/>
      <value value="0.4"/>
      <value value="0.5"/>
      <value value="0.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="2.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="added-d-juv-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110406c Calibrate k and e-use" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="m-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.3"/>
      <value value="0.4"/>
      <value value="0.5"/>
      <value value="0.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="2.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="added-d-juv-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110406d Calibrate k and e-use" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="m-survival-prob-const">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.3"/>
      <value value="0.4"/>
      <value value="0.5"/>
      <value value="0.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="2.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="added-d-juv-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110406e Calibrate k and e-use" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="m-survival-prob-const">
      <value value="0.25"/>
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="2.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="added-d-juv-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110406f Calibrate k and e-use" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="m-survival-prob-const">
      <value value="0.5"/>
      <value value="0.75"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="2.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="added-d-juv-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110411 Calibrate k and e-use" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="m-survival-prob-const">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="2.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="added-d-juv-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110412a Calibrate juv m" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-juv-mort-const">
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="2.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110412b Calibrate juv m" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-juv-mort-const">
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="2.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110412c Calibrate juv m" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-juv-mort-const">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="2.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110412d Calibrate e-use" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="2.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-juv-mort-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110412e Calibrate e-use (REFERENCE)" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-juv-mort-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110412f Calibrate e-use" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-juv-mort-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110415 - Reference model" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-juv-mort-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110415 - Reference + wind farms" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;All&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-juv-mort-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110415 - Reference + ships" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-juv-mort-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110415 - Reference + ships + windf" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;All&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-juv-mort-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110415 - Reference + pot wind farms" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;User-def&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-juv-mort-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110415 - Reference + pot wf constr" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;User-def&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-juv-mort-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="true"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110415 - Reference model + extra juv mort" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-juv-mort-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="0.9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110506 - Reference model + bycatch-025" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-juv-mort-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110506&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.025"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110506 - Reference model + bycatch-050" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-juv-mort-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110506&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110506 - Reference model + bycatch-100" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-juv-mort-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110506&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110506 - Reference model + bycatch-000" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-juv-mort-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110506&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="111007 - Calibrate adult mort m" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-juv-mort-const">
      <value value="0.5"/>
    </enumeratedValueSet>
    <steppedValueSet variable="m-mort-prob-const" first="0.3" step="0.1" last="0.7"/>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="1.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-110404&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="111129 calibrate x-surv-const and e-use-1" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.6"/>
      <value value="0.1"/>
      <value value="0.9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="5.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111126&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="111129 calibrate x-surv-const and e-use-2" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.6"/>
      <value value="0.1"/>
      <value value="0.9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111126&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="111129 calibrate x-surv-const and e-use-3" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 0</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.6"/>
      <value value="0.1"/>
      <value value="0.9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="6.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111126&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="111129 calibrate x-surv-const and e-use-4" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.6"/>
      <value value="0.1"/>
      <value value="0.9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111126&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="111129 calibrate x-surv-const and e-use-5" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.6"/>
      <value value="0.1"/>
      <value value="0.9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111126&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="111129 calibrate x-surv-const and e-use-6" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.6"/>
      <value value="0.1"/>
      <value value="0.9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111126&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="111129 calibrate x-surv-const and e-use-7" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.6"/>
      <value value="0.1"/>
      <value value="0.9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="3.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111126&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="111129 calibrate x-surv-const and e-use-8" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111126&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="111129 calibrate x-surv-const and e-use-9" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111126&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="111220 New reference, no bycatch" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="111220 New reference + all wind farms" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;All&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="111129 calibrate x-surv-const and e-use-10" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="3.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111126&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="111129 calibrate x-surv-const and e-use-11" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111126&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="111220 New reference + all wind farms + all ships" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;All&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="111220 New reference + all wind farms + pot wf" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;User-def&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="111220 New reference + all wind farms + pot wf + all ships" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;User-def&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120124 New reference + all wind farms + bycatch=0.01" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;All&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120124 New reference + all wind farms + bycatch=0.017" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;All&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120124 New reference + all wind farms + bycatch=0.05" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;All&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120124 New reference + all wind farms + bycatch=0.5" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;All&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120309 New reference, no bycatch, but rU=0.02" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120309 New reference + all wind farms + pot wf + all ships but rU=0.02" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;User-def&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120309 New reference + all wind farms + pot wf but rU=0.02" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;User-def&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120309 New reference + all wind farms + bycatch=0.02" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;All&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120323 New reference, no bycatch, but rU=0.06" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.06"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120323 New reference + all wind farms + pot wf + all ships but rU=0.06" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.06"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;User-def&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120323 New reference + all wind farms + pot wf but rU=0.06" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.06"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;User-def&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120323 New reference + all wind farms + bycatch=0.02" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;All&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120402 New reference, no bycatch, but rU=0.08" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.08"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120402 New reference + all wind farms + pot wf + all ships but rU=0.08" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.08"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;User-def&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120402 New reference + all wind farms + pot wf but rU=0.08" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.08"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;User-def&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120402 New reference, no bycatch, but rU=0.09" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.09"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120402 New reference + all wind farms + pot wf + all ships but rU=0.09" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.09"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;User-def&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120402 New reference + all wind farms + pot wf but rU=0.09" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.09"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;User-def&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120419 New reference + all wind farms + bycatch=0.1" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;All&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120419 New reference + all wind farms + bycatch=0.2" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;All&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120503 New reference, no bycatch, but rU=0.005" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.0050"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120508 New reference + all wind farms + pot wf + all ships but rU=0.08" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.08"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;User-def&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120508 New reference + all wind farms but rU=0.09" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.09"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;All&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120508 New reference + all wind farms + all ships but rU=0.09" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.09"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;All&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120521 New reference + pot wf + all ships" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;Potential&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120521 New reference + bycatch=0.02" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120521 New reference + bycatch=0.05" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120612 New reference + all ships" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120612 New reference + all wind farms + all ships" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;All&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120612 New reference + potential wind farms" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deter-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;Potential&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120708 New reference + all wind farms + pot wf, dt=4" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="deter-time">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;User-def&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120708 New reference + all wind farms + pot wf, dt=3" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="deter-time">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;User-def&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120708 New reference + all wind farms + pot wf, dt=2" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="deter-time">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;User-def&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120708 New reference + all wind farms + pot wf, dt=1" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-output word random 10000 ".txt"</final>
    <timeLimit steps="691200"/>
    <exitCondition>count porps = 1</exitCondition>
    <metric>count porps</metric>
    <metric>mean [energy-level] of porps</metric>
    <metric>sum [ food-level ] of patches with [ food-level &gt; 0 ]</metric>
    <enumeratedValueSet variable="deter-time">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="x-survival-prob-const">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-30min">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim-day">
      <value value="14400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-mort-prob-const">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deterrence-coeff">
      <value value="7.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="&quot;monthly&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-deterrence-dist">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disp-depth">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-disp-targets">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pile-driving">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bycatch-prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-dist-to-target">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output-name">
      <value value="&quot;kat-111220&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;Kattegat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="wind-farms">
      <value value="&quot;User-def&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-disp-dist">
      <value value="1.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incl-ships">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="e-use-per-km">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180

@#$#@#$#@
0
@#$#@#$#@
