; Changes made 6/19/24
;;Built from Pulmonary_Histology_Python_7302023_colexperiments
;;Edited from Pulmonary_Histology_5_20_2024_Bill
;; This version has more specific collagen deposition amounts

extensions [ py ]

breed [ thy1ps thy1p ] ; Thy1 positive - Orange
breed [ thy1ns thy1n ] ; Thy1 negative - Blue
breed [ gen-fibs gen-fib ] ; Generic Fibroblasts prior to activation - Red

globals [
  gen-fib-count ; Total number of generic fibroblasts in the beginning
  p-collagen-deposition ; Rate of collagen dep. for thy1p
  n-collagen-deposition ; Rate of collagen dep. for thy1n
  collagen-degradation ; Rate of collagen degradation
 ; AS-entry-threshold ; Threshold to enter alveolar space
  AS-entry-random ; Stochastic factor to enter alveolar space

  daily-migration-rate
  ave-lifespan
  dev-lifespan
  prolif-rate
  apop-rate

  tnf-thresh
  il1b-thresh

  colcontent
  fibrosis-score

  taval
  pcolplus

  alv-count
  final-thy1p
  final-thy1n

  image-title
  image-mask-title

  mask-psource
  mask-nsource
  mask-type

]

patches-own [
 matrix
 tissue-type
 ;pcolplus
 c
 ptnfa
 pil1b
 tgfbenv
 rand
 pat-psource
 pat-nsource
 pat-mask
]

turtles-own [
 age
 convcount
  tgfbin
  mechin
  thy1in
  fibronin
  faslin
  fiboutput
  asma
  col1
  col3
  prolif
  apop
  colplus
  tnfa
  il1b
  psource
  nsource
  myid
  rancol
]

to setup
  clear-all
  reset-ticks

  setup-globals
  setup-patches
  setup-thy1ps
  setup-thy1ns
  set TGFBint tgfbint
;    py:setup py:python
;  (py:run
;    "import numpy as np"
;    "import os"
;    "from fibroblast_v3_TCBoth_NetfluxODE_run_NetLogo import oderun"
;    )
end

to setup-globals
  set gen-fib-count thy1ps-count + thy1ns-count
  set p-collagen-deposition 0.05
  set n-collagen-deposition 0.25
  set collagen-degradation 0.5
  set AS-entry-threshold 4
  set AS-entry-random 40
  set daily-migration-rate 1
  set ave-lifespan 57
  set dev-lifespan 3
  set prolif-rate 0.018
  set apop-rate 0.072
  set tnf-thresh set-tnf-thresh
  set il1b-thresh set-il1b-thresh
  set mask-psource 0
  set mask-nsource 0

end


to setup-patches

  ; Import the picture of pulmonary tissue

  ask patches [ set pcolor white ]
  import-pcolors "healthylung_blue.jpg"
  ask patches [
  set matrix  0.1 ;matrix works as a counter to detemine patch stiffness
  ]

;  ; Create contrast from the tissue data
  ask patches [
    ifelse pcolor < 90 [
      set pcolor white
      set tissue-type "alveolar-space"
    ] [
      set pcolor grey
      set tissue-type "interstitial"
    ]
      set c count neighbors with [pcolor = white]
      if c >= 4 [set pcolor white set tissue-type "alveolar-space"]
  ]
  ask patches [
      set c count neighbors with [pcolor = white]
      if c >= 8 [set pcolor white set tissue-type "alveolar-space"]
  ]
  ask patches with [ tissue-type = "interstitial" ]  [
    set pil1b 50
    set ptnfa 50
  ]

    ask patches [
    set rand ((random-float 2.00000001))
    if rand < 0.2 [set rand 0.2]
    set rand precision rand 3
    set tgfbenv (tgfbint * rand)

    set pat-psource 0
    set pat-nsource 0
    set pcolplus 0
  ]


end


; DEPRECATED: Use this to start the board with a random number of thy1ps in interstitial space
to setup-thy1ps

  let num-interstitial count patches with [ tissue-type = "interstitial" ]

  while [ count thy1ps < thy1ps-count ] [
    ask patches with [ tissue-type = "interstitial" ]  [
      if count thy1ps < thy1ps-count [
        if random num-interstitial < thy1ps-count [ sprout-thy1ps 1 ]
      ]
    ]
  ]

  ask thy1ps [ set shape "cell" ]
  ask thy1ps [ set color 26]
  ask thy1ps [ set size 3]

  ask thy1ps [ set age  random 25 ]
  ask thy1ps [ set convcount 0 ]

  ask thy1ps [
    set mechin 0.1
    set tgfbin 0.05
    set thy1in 1
    set fibronin 1
    set faslin 0.25
    set tnfa random 10
    set il1b random 10
    set asma 0.0015
    set col1 0.028398702
    set col3 0.028398702
    set prolif 0.0003
    set apop 0.0008
    set psource who
    set nsource 0
    set myid who
    set label "" ;psource
    set label-color black
  ]



end

; DEPRECATED: Use this to start the board with a random number of thy1ns in interstitial space
to setup-thy1ns

  let num-interstitial count patches with [ tissue-type = "interstitial" ]

  while [ count thy1ns < thy1ns-count ] [
    ask patches with [ tissue-type = "interstitial" ]  [
      if count thy1ns < thy1ns-count [
        if random num-interstitial < thy1ns-count [ sprout-thy1ns 1 ]
      ]
    ]
  ]


  ask thy1ns [ set shape "cell" ]
  ask thy1ns [ set color 85]
  ask thy1ns [ set size 3]

  ask thy1ns [ set age random 25 ]
  ask thy1ns [ set convcount 0 ]

    ask thy1ns [
    set mechin 0.1
    set tgfbin 0.05
    set thy1in 0.25
    set fibronin 1
    set faslin 0.25
    set tnfa 0
    set il1b 0
    set asma 0.0012
    set col1 0.028205855
    set col3 0.028205855
    set prolif 0.0003
    set apop 0
    set psource 0
    set nsource who
    set myid who
    set label "" ;nsource
    set label-color yellow
  ]


end

;;; have mask that are agents of different families linked and if its linked to a negative the link is a different color


to go
  ;if ticks >= end-tick [
   ;print "end"
   calc-fib-score
   ;set image-title (word "results "  random-float 1.0 ".png")
   ;export-view (image-title) ;; Colon characters in the time cause errors on Windows
   ;hide-turtles
   ;create-mask
   ;set image-mask-title (word "mask- " image-title ".png")
   ;export-view image-mask-title

  ; stop
 ; ]  ;; stop after 730 ticks
  move-thy1ps
  move-thy1ns
  uptake_tgfb

 ; runpy
  fibcheck
  colldepon

  stiffen

  age-fib
  prolif-fib
  convert-fibs

  tick                    ;; increase the tick counter by 1 each time through
end

to runpy
  ask turtles [
    if breed != gen-fibs [
    ;tgfbin,mechin,il6in,il1in,tnfain,pdgfin,fgfin,hypoxia,thy1state
  py:set "tgfbin" tgfbin
  py:set "mechin" mechin
  py:set "thy1in" thy1in
  py:set "fibronin" fibronin
  py:set "faslin" faslin
  set fiboutput (
     py:runresult "oderun(mechin,tgfbin,thy1in,fibronin,faslin)"
     )
   set asma item 0 fiboutput
   set col1 item 1 fiboutput
   set col3 item 2 fiboutput
   set prolif item 3 fiboutput
   set apop item 4 fiboutput
    ]
  ]

end

; Increase the age with every tick
to age-fib
  ask turtles [ set age age + 1 ]
  ask turtles [ set convcount convcount + 1 ]
;  ask turtles [
;    if age >= ave-lifespan + random dev-lifespan[
;      ask thy1ps [
;    if ((apop *  random 100)) > apop-rate [
;          die]
;      ]
;      ask thy1ns [
;            if (((apop + random-float 0.001) * random 100)) > apop-rate [
;          die]
;      ]
;  ]
;  ]

  ask turtles [
    if ((apop * 1000) * ( random dev-lifespan / 2)) + (ave-lifespan - (dev-lifespan / 2)) < age [ die ] ; If they are within the average lifespan range, with a stochastic factor, they die
  ]

end

; Randomly proliferate at some point in their life
to prolif-fib
  ask turtles [
   if count turtles < 150 [
    if (prolif * random 100)  > prolif-rate [
        hatch 1 [
          set age 0
          set convcount 0
          set mechin 0.25
          set tgfbin 0.05
          set tnfa 0
          set il1b 0
          set myid who
          ifelse nsource = 0
          [set label "" ;psource
          set label-color black]
          [set label "" ;nsource
          set label-color yellow]
        ]
        ask one-of turtles [
          if age >= ave-lifespan + dev-lifespan [
          die]
        ]

  ]
  ]
  ]
end

to move-thy1ps ; agent turns to a random angle then moves forward one patch
  ask thy1ps [
    right random 360
    if patch-ahead 1 != nobody [
    if [pcolor] of patch-ahead 1 != white
      [ forward daily-migration-rate ]
      if [pcolor] of patch-here = 125
      [ forward daily-migration-rate]
    set il1b il1b + random 5
    set tnfa tnfa + random 5
    ]

  ]
end

to move-thy1ns ; agent turns to a random angle then moves forward one patch

  ask thy1ns [
    right random 360
    if patch-ahead 1 != nobody [
    ; Only thy1 negative fibs can enter the alveolar space
    if [tissue-type] of patch-ahead 1 = "interstitial" and [tissue-type] of patch-ahead 2 = "interstitial" and [pcolor] of patch-ahead 1 != white and [pcolor] of patch-ahead 2 != white [ forward daily-migration-rate ]
    if [tissue-type] of patch-ahead 1 = "alveolar-space" [
      if [matrix] of patch-at 1 0 + [matrix] of patch-at -1 0 + [matrix] of patch-at 0 1 + [matrix] of patch-at 0 -1 + [matrix] of patch-here >= AS-entry-threshold  [
        ;print("Enter alveolar space. ")
        forward daily-migration-rate
      ]
      if [pcolor] of patch-here = 125
      [ forward daily-migration-rate]
    ]
    ]
  ]
end

to fibcheck       ;Thy-1 positive cells check for collagen on patch ahead
  ask thy1ps [
    if patch-ahead 1 != nobody[
;    if [pcolor] of patch-here != grey [
;      set col1 0.0617
;      set col3 0.0617
;    ]
      if [matrix] of patch-here < 0.29 [
        set col1 0.028398702
        set col3 0.028398702] ;0.2
      if [matrix] of patch-here > 0.3 [
        set col1 0.040149463
        set col3 0.040149463] ;0.3
      if [matrix] of patch-here > 0.4 [
        set col1 0.056530381
        set col3 0.056530381] ;0.4
      if [matrix] of patch-here > 0.5 [
        set col1 0.077850997
        set col3 0.077850997] ;0.5
      if [matrix] of patch-here > 0.6 [
        set col1 0.104572381
        set col3 0.104572381] ;0.6
      if [matrix] of patch-here > 0.7 [
        set col1 0.137367574
        set col3 0.137367574] ;0.7
      if [matrix] of patch-here > 0.8 [
        set col1 0.177238071
        set col3 0.177238071] ;0.8
      if [matrix] of patch-here > 0.9 [
        set col1 0.225733939
        set col3 0.225733939] ;0.9

    set colplus (((col1 + col3) * coldep) - coldeg * abs(matrix))
    set pcolplus colplus
    set mechin matrix
    if mechin > 1 [
     set mechin 1]

    ifelse [pcolor] of patch-ahead 1 != white and [pcolor] of patch-ahead 1 != grey [          ;First they check for "stiff" patches, if patch is stiff lay down more matrix

      ask patch-ahead 1 [
   set matrix (matrix + pcolplus)
          if matrix > 1 [set matrix 1]
      ]
      ]
    [set rancol (0 + random 1000)
        if rancol = 1000 [
        set colplus (((col1 + col3) * coldep) - coldeg * matrix)
        set pcolplus colplus
        set mechin matrix
        if mechin > 1 [
          set mechin 1
        ]
      ask patch-ahead 1 [
   set matrix (matrix + pcolplus)
            if matrix > 1 [set matrix 1]
        ]
        ]
      ]
  ]
  ]
end

to colldepon
  ask thy1ns [             ;Thy-1 negative fibroblasts lay down collagen
    if patch-ahead 1 != nobody[
;      if [pcolor] of patch-here != grey and [pcolor] of patch-here != white [
;      set col1 0.0365
;      set col3 0.0365
;      ]
      if [matrix] of patch-here < 0.29 [
        set col1 0.028205855
        set col3 0.028205855] ;0.2
      if [matrix] of patch-here > 0.3 [
        set col1 0.032607399
        set col3 0.032607399] ;0.3
      if [matrix] of patch-here > 0.4 [
        set col1 0.0386445
        set col3 0.0386445] ;0.4
      if [matrix] of patch-here > 0.5 [
        set col1 0.046444017
        set col3 0.046444017] ;0.5
      if [matrix] of patch-here > 0.6 [
        set col1 0.056133012
        set col3 0.056133012] ;0.6
      if [matrix] of patch-here > 0.7 [
        set col1 0.067870448
        set col3 0.067870448] ;0.7
      if [matrix] of patch-here > 0.8 [
        set col1 0.081930506
        set col3 0.081930506] ;0.8
      if [matrix] of patch-here > 0.9 [
        set col1 0.098841275
        set col3 0.098841275] ;0.9
    set colplus (((col1 + col3) * coldep)- coldeg * abs(matrix) )
    set pcolplus colplus
    set mechin matrix
      if mechin > 1 [
      set mechin 1
      ]
  ask patch-ahead 1 [
    set matrix (matrix + pcolplus)
        if matrix > 1 [set matrix 1]
  ]
  ]
  ]
end

to stiffen   ;patches change color and "stiffness" based on matrix value
  ask patches [
    if matrix > 1 [ set matrix 1 ]
    if matrix < 0 [set matrix 0]
    if matrix > 0.2 [ set pcolor 138]
    if matrix > 0.3 [set pcolor 137]
    if matrix > 0.4 [set pcolor 136]
    if matrix > 0.5 [set pcolor 135
     set ptnfa ptnfa + 50
     set ptnfa ptnfa + 50]
    if matrix > 0.6 [set pcolor 128]
    if matrix > 0.7 [set pcolor 127
        set ptnfa ptnfa + 50
        set ptnfa ptnfa + 50
    ]
    if matrix > 0.8 [set pcolor 126]
    if matrix > 0.9 [set pcolor 125]

;    if pcolor = 128 [   ;Become "stiff"
;    if matrix > 0.75
;    [set pcolor 125
;     set ptnfa ptnfa + 50
;     set ptnfa ptnfa + 50
;
;      ]
;    ]
;
;    if pcolor = grey or pcolor = white [    ;Become "medium"
;    if matrix > 0.5
;       [set pcolor 128
;        set ptnfa ptnfa + 50
;        set ptnfa ptnfa + 50
;
;      ]
;    ]
    ]
end


;Thy1 positive fibroblasts (thy1ps) convert to thy1n fibroblasts over time
to convert-fibs

  ask thy1ps [
;    if transition? = true [
      if (count thy1ns / count turtles) < transition-threshold [
    if tnfa > tnf-thresh and il1b > il1b-thresh [
    set breed thy1ns
    set shape "cell"
    set color 85
    set thy1in 0.25
    set asma 0.0012
    set col1 0.001
    set col3 0.001
    set prolif 0.0003
    set apop 0
      ]
    ]
  ]
;  ]
end

to uptake_tgfb
  ask turtles [
  rt random 360
  if patch-ahead 1 != nobody [
  ;tgfb uptake
    if [tgfbenv] of patch-ahead 1 > 2
      [set tgfbin (tgfbenv / ( tgfbenv + 700 ) )
        set taval tgfbin
         ask patch-ahead 1 [
          set tgfbenv (tgfbenv - (taval * random 0.05))]
        if tgfbin > 1 [ set tgfbin 1]
  ]

    ]
  ]
end

to calc-fib-score
  set colcontent ((count patches with [pcolor != white and pcolor != grey]) / count patches) * 100
  if colcontent < 25
  [set fibrosis-score 1]
  if colcontent > 25 and colcontent < 75
  [set fibrosis-score 2]
  if colcontent > 75
  [set fibrosis-score 3]
  ;print fibrosis-score

  set alv-count count patches with [tissue-type = "alveolar-space" and pcolor != white] / count patches with [tissue-type = "alveolar-space"]
  set final-thy1p count thy1ps
  set final-thy1n count thy1ns
end

to create-mask
  ask patches [set pcolor white]
  ask turtles [
        set mask-psource psource
        set mask-nsource nsource
    if thy1in = 1 [set mask-type 1]
    if thy1in != 1 [set mask-type 0]
    ;hide-turtle
    ask patch-here
    [set pat-psource mask-psource
     set pat-nsource mask-nsource
     set pat-mask mask-type
    ]
    ;hide-turtle
  ]
  ask patches   [
    if pat-psource != 0 [
    set pcolor green
    set plabel pat-psource
    if pat-mask = 1
      [set plabel-color black]
    if pat-mask = 0
      [set plabel-color orange]
    ]
    if pat-nsource != 0 [
    set pcolor pink
      ;set pcolor scale-color red pat-nsource 50 0
    set plabel pat-nsource
    set plabel-color black]
  ]
end
@#$#@#$#@
GRAPHICS-WINDOW
703
22
1286
486
-1
-1
5.0
1
10
1
1
1
0
1
1
1
-57
57
-45
45
1
1
1
ticks
30.0

SLIDER
13
54
185
87
thy1ps-count
thy1ps-count
0
100
90.0
1
1
NIL
HORIZONTAL

SLIDER
13
109
185
142
thy1ns-count
thy1ns-count
0
100
10.0
1
1
NIL
HORIZONTAL

BUTTON
17
15
81
48
Setup
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
92
15
155
48
Go
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

MONITOR
217
53
288
98
Thy 1 Pos
count thy1ps
17
1
11

MONITOR
217
107
287
152
Thy 1 Neg
count thy1ns
17
1
11

MONITOR
301
55
384
100
Soft Patches
count patches with [pcolor = white or pcolor = 95]
17
1
11

MONITOR
303
109
405
154
Medium Patches
count patches with [pcolor = 128]
17
1
11

MONITOR
303
166
386
211
Stiff Patches
count patches with [pcolor = 125]
17
1
11

PLOT
12
218
307
399
Fibroblast Count over Time
Time
Number of Fibroblasts
0.0
100.0
0.0
10.0
true
true
"" ""
PENS
"Thy 1 positive" 1.0 0 -14439633 true "" "plot count thy1ps"
"Thy 1 negative" 1.0 0 -16777216 true "" "plot count thy1ns"

PLOT
413
53
689
203
Surface Stiffness over Time
Time
Different Stiffness Patches
0.0
100.0
0.0
1100.0
true
true
"" ""
PENS
"Soft " 1.0 0 -13840069 true "" "plot count patches with [pcolor = white]"
"Medium" 1.0 0 -987046 true "" "plot count patches with [pcolor = 128]"
"Stiff" 1.0 0 -2674135 true "" "plot count patches with [pcolor = 125]"

MONITOR
557
489
668
534
% AV w Collagen 
count patches with [tissue-type = \"alveolar-space\" and pcolor != white] / count patches with [tissue-type = \"alveolar-space\"]
17
1
11

MONITOR
445
490
544
535
AV w Collagen
count patches with [tissue-type = \"alveolar-space\" and pcolor != white]
17
1
11

PLOT
319
216
617
398
Alveolar Patches with Collagen over Time 
Time 
Patches 
0.0
1000.0
0.0
1000.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot count patches with [tissue-type = \"alveolar-space\" and pcolor != white]"

SLIDER
64
413
236
446
set-tnf-thresh
set-tnf-thresh
0
100
50.0
1
1
NIL
HORIZONTAL

SLIDER
290
412
462
445
set-il1b-thresh
set-il1b-thresh
0
100
50.0
1
1
NIL
HORIZONTAL

SLIDER
281
449
481
482
transition-threshold
transition-threshold
0
1
0.1
.1
1
NIL
HORIZONTAL

SLIDER
63
491
235
524
coldep
coldep
0
10000
100.0
.001
1
NIL
HORIZONTAL

MONITOR
192
165
296
210
NIL
fibrosis-score
0
1
11

SLIDER
252
494
424
527
AS-entry-threshold
AS-entry-threshold
0
10
4.0
.1
1
NIL
HORIZONTAL

SLIDER
62
535
234
568
coldeg
coldeg
0
200
5.0
0.001
1
NIL
HORIZONTAL

SLIDER
65
448
237
481
TGFBint
TGFBint
0
1400
700.0
1
1
pg/mL
HORIZONTAL

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
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

cell
true
13
Polygon -2064490 true true 135 0 15 90 60 270 240 270 285 210 270 75 180 60
Circle -5825686 true false 105 105 90

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

orbit 1
true
0
Circle -7500403 true true 116 11 67
Circle -7500403 false true 41 41 218

orbit 2
true
0
Circle -7500403 true true 116 221 67
Circle -7500403 true true 116 11 67
Circle -7500403 false true 44 44 212

orbit 3
true
0
Circle -7500403 true true 116 11 67
Circle -7500403 true true 26 176 67
Circle -7500403 true true 206 176 67
Circle -7500403 false true 45 45 210

orbit 4
true
0
Circle -7500403 true true 116 11 67
Circle -7500403 true true 116 221 67
Circle -7500403 true true 221 116 67
Circle -7500403 false true 45 45 210
Circle -7500403 true true 11 116 67

orbit 5
true
0
Circle -7500403 true true 116 11 67
Circle -7500403 true true 13 89 67
Circle -7500403 true true 178 206 67
Circle -7500403 true true 53 204 67
Circle -7500403 true true 220 91 67
Circle -7500403 false true 45 45 210

orbit 6
true
0
Circle -7500403 true true 116 11 67
Circle -7500403 true true 26 176 67
Circle -7500403 true true 206 176 67
Circle -7500403 false true 45 45 210
Circle -7500403 true true 26 58 67
Circle -7500403 true true 206 58 67
Circle -7500403 true true 116 221 67

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

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

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

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.2.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="Vary-Thy1p" repetitions="5" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="1000"/>
    <metric>count patches with [pcolor = 128] / count patches</metric>
    <metric>count patches with [pcolor = 125] / count patches</metric>
    <metric>count patches with [tissue-type = "alveolar-space" and pcolor != white] / count patches with [tissue-type = "alveolar-space"]</metric>
    <enumeratedValueSet variable="collagen-degradation">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="degrader-count">
      <value value="10"/>
    </enumeratedValueSet>
    <steppedValueSet variable="thy1ps-count" first="0" step="10" last="100"/>
    <enumeratedValueSet variable="n-collagen-deposition">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AS-entry-threshold">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="thy1ns-count">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-collagen-deposition">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AS-entry-random">
      <value value="50"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Vary-Thy1n" repetitions="4" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="730"/>
    <metric>count patches with [pcolor = 128] / count patches</metric>
    <metric>count patches with [pcolor = 125] / count patches</metric>
    <metric>count patches with [tissue-type = "alveolar-space" and pcolor != white] / count patches with [tissue-type = "alveolar-space"]</metric>
    <enumeratedValueSet variable="thy1ps-count">
      <value value="0"/>
      <value value="15"/>
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="thy1ns-count">
      <value value="0"/>
      <value value="15"/>
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-conv">
      <value value="365"/>
      <value value="800"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Vary-Degraders" repetitions="5" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="1000"/>
    <metric>count patches with [pcolor = 128] / count patches</metric>
    <metric>count patches with [pcolor = 125] / count patches</metric>
    <metric>count patches with [tissue-type = "alveolar-space" and pcolor != white] / count patches with [tissue-type = "alveolar-space"]</metric>
    <enumeratedValueSet variable="collagen-degradation">
      <value value="0.5"/>
    </enumeratedValueSet>
    <steppedValueSet variable="degrader-count" first="0" step="10" last="100"/>
    <enumeratedValueSet variable="thy1ps-count">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-collagen-deposition">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AS-entry-threshold">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="thy1ns-count">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-collagen-deposition">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AS-entry-random">
      <value value="50"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="AS-Entry-Bivariate" repetitions="5" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="1000"/>
    <metric>count patches with [pcolor = 128] / count patches</metric>
    <metric>count patches with [pcolor = 125] / count patches</metric>
    <metric>count patches with [tissue-type = "alveolar-space" and pcolor != white] / count patches with [tissue-type = "alveolar-space"]</metric>
    <enumeratedValueSet variable="collagen-degradation">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="degrader-count">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="thy1ps-count">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-collagen-deposition">
      <value value="1"/>
    </enumeratedValueSet>
    <steppedValueSet variable="AS-entry-threshold" first="0" step="3" last="30"/>
    <enumeratedValueSet variable="thy1ns-count">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-collagen-deposition">
      <value value="0.5"/>
    </enumeratedValueSet>
    <steppedValueSet variable="AS-entry-random" first="0" step="10" last="100"/>
  </experiment>
  <experiment name="Global-Sampling" repetitions="3" sequentialRunOrder="false" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="1000"/>
    <metric>count patches with [pcolor = 128] / count patches</metric>
    <metric>count patches with [pcolor = 125] / count patches</metric>
    <metric>count patches with [tissue-type = "alveolar-space" and pcolor != white] / count patches with [tissue-type = "alveolar-space"]</metric>
    <enumeratedValueSet variable="collagen-degradation">
      <value value="0.5"/>
    </enumeratedValueSet>
    <steppedValueSet variable="degrader-count" first="0" step="10" last="30"/>
    <steppedValueSet variable="thy1ps-count" first="0" step="10" last="30"/>
    <enumeratedValueSet variable="n-collagen-deposition">
      <value value="1"/>
    </enumeratedValueSet>
    <steppedValueSet variable="AS-entry-threshold" first="0" step="6" last="30"/>
    <steppedValueSet variable="thy1ns-count" first="0" step="10" last="30"/>
    <enumeratedValueSet variable="p-collagen-deposition">
      <value value="0.5"/>
    </enumeratedValueSet>
    <steppedValueSet variable="AS-entry-random" first="0" step="25" last="100"/>
  </experiment>
  <experiment name="Convertion Experiment" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="730"/>
    <metric>count thy1ps</metric>
    <metric>count thy1ns</metric>
    <metric>count degraders</metric>
    <metric>count patches with [pcolor = 128] / count patches</metric>
    <metric>count patches with [pcolor = 125] / count patches</metric>
    <metric>count patches with [tissue-type = "alveolar-space" and pcolor != white] / count patches with [tissue-type = "alveolar-space"]</metric>
    <enumeratedValueSet variable="degrader-count">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="thy1ps-count">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="thy1ns-count">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="collagen-degradation">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-collagen-deposition">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AS-entry-threshold">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-collagen-deposition">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AS-entry-random">
      <value value="50"/>
    </enumeratedValueSet>
    <steppedValueSet variable="set-conv" first="0" step="100" last="800"/>
    <enumeratedValueSet variable="set-conv2">
      <value value="500"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="maxxed-Vary-Thy1n" repetitions="20" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="730"/>
    <metric>count patches with [pcolor = 128]</metric>
    <metric>count patches with [pcolor = 125]</metric>
    <metric>count patches with [tissue-type = "alveolar-space" and pcolor != white] / count patches with [tissue-type = "alveolar-space"]</metric>
    <metric>count patches</metric>
    <enumeratedValueSet variable="degrader-count">
      <value value="0"/>
      <value value="6"/>
      <value value="9"/>
      <value value="12"/>
      <value value="18"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="thy1ps-count">
      <value value="0"/>
      <value value="6"/>
      <value value="9"/>
      <value value="12"/>
      <value value="18"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="thy1ns-count">
      <value value="0"/>
      <value value="12"/>
      <value value="18"/>
      <value value="24"/>
      <value value="36"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="collagen-degradation">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-collagen-deposition">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AS-entry-threshold">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-collagen-deposition">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AS-entry-random">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-conv">
      <value value="183"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Population Exp" repetitions="10" sequentialRunOrder="false" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="730"/>
    <metric>count patches with [pcolor = 128] / count patches</metric>
    <metric>count patches with [pcolor = 125] / count patches</metric>
    <metric>count patches with [tissue-type = "alveolar-space" and pcolor != white] / count patches with [tissue-type = "alveolar-space"]</metric>
    <enumeratedValueSet variable="thy1ps-count">
      <value value="0"/>
      <value value="15"/>
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="thy1ns-count">
      <value value="0"/>
      <value value="15"/>
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-conv">
      <value value="55"/>
      <value value="356"/>
      <value value="500"/>
      <value value="700"/>
      <value value="800"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Vary-All" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="1095"/>
    <metric>count patches with [pcolor = 128] / count patches</metric>
    <metric>count patches with [pcolor = 125] / count patches</metric>
    <metric>count patches with [tissue-type = "alveolar-space" and pcolor != white] / count patches with [tissue-type = "alveolar-space"]</metric>
    <enumeratedValueSet variable="thy1ps-count">
      <value value="0"/>
      <value value="25"/>
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="thy1ns-count">
      <value value="0"/>
      <value value="25"/>
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-tnf-thresh">
      <value value="30"/>
      <value value="60"/>
      <value value="90"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-il1b-thresh">
      <value value="30"/>
      <value value="60"/>
      <value value="90"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Vary-Populations Transitioning" repetitions="5" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="1095"/>
    <metric>count patches with [pcolor = 128] / count patches</metric>
    <metric>count patches with [pcolor = 125] / count patches</metric>
    <metric>count patches with [tissue-type = "alveolar-space" and pcolor != white] / count patches with [tissue-type = "alveolar-space"]</metric>
    <enumeratedValueSet variable="thy1ps-count">
      <value value="0"/>
      <value value="5"/>
      <value value="25"/>
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="thy1ns-count">
      <value value="0"/>
      <value value="5"/>
      <value value="25"/>
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-tnf-thresh">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-il1b-thresh">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="transition?">
      <value value="false"/>
      <value value="true"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Vary-Populations Transitioning" repetitions="5" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="1095"/>
    <metric>count patches with [pcolor = 128] / count patches</metric>
    <metric>count patches with [pcolor = 125] / count patches</metric>
    <metric>count patches with [tissue-type = "alveolar-space" and pcolor != white] / count patches with [tissue-type = "alveolar-space"]</metric>
    <enumeratedValueSet variable="thy1ps-count">
      <value value="25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="thy1ns-count">
      <value value="25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-tnf-thresh">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-il1b-thresh">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="transition?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="transition-threshold">
      <value value="0.4"/>
      <value value="0.8"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="10percent Low Ceiling Transitioning" repetitions="4" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="1095"/>
    <metric>count patches with [pcolor = 128] / count patches</metric>
    <metric>count patches with [pcolor = 125] / count patches</metric>
    <metric>count patches with [tissue-type = "alveolar-space" and pcolor != white] / count patches with [tissue-type = "alveolar-space"]</metric>
    <metric>count thy1ps</metric>
    <metric>count thy1ns</metric>
    <enumeratedValueSet variable="thy1ps-count">
      <value value="0"/>
      <value value="5"/>
      <value value="45"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="thy1ns-count">
      <value value="0"/>
      <value value="5"/>
      <value value="45"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-tnf-thresh">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-il1b-thresh">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="transition?">
      <value value="true"/>
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="transition-threshold">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="10percent No Transitioning" repetitions="3" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="1095"/>
    <metric>count patches with [pcolor = 128] / count patches</metric>
    <metric>count patches with [pcolor = 125] / count patches</metric>
    <metric>count patches with [tissue-type = "alveolar-space" and pcolor != white] / count patches with [tissue-type = "alveolar-space"]</metric>
    <metric>count thy1ps</metric>
    <metric>count thy1ns</metric>
    <enumeratedValueSet variable="thy1ps-count">
      <value value="45"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="thy1ns-count">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-tnf-thresh">
      <value value="30"/>
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-il1b-thresh">
      <value value="30"/>
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="transition?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="transition-threshold">
      <value value="0.1"/>
      <value value="0.2"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="1yr_transition_thresh50_multicase" repetitions="5" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="365"/>
    <metric>count patches with [pcolor = 128] / count patches</metric>
    <metric>count patches with [pcolor = 125] / count patches</metric>
    <metric>count patches with [tissue-type = "alveolar-space" and pcolor != white] / count patches with [tissue-type = "alveolar-space"]</metric>
    <metric>count thy1ps</metric>
    <metric>count thy1ns</metric>
    <enumeratedValueSet variable="thy1ps-count">
      <value value="45"/>
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="thy1ns-count">
      <value value="0"/>
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-tnf-thresh">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-il1b-thresh">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="transition?">
      <value value="true"/>
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="transition-threshold">
      <value value="0.1"/>
      <value value="0.2"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Colvar_Testing_v2" repetitions="5" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>count patches with [pcolor = 128] / count patches</metric>
    <metric>count patches with [pcolor = 125] / count patches</metric>
    <metric>count patches with [tissue-type = "alveolar-space" and pcolor != white] / count patches with [tissue-type = "alveolar-space"]</metric>
    <metric>count thy1ps</metric>
    <metric>count thy1ns</metric>
    <metric>fibrosis-score</metric>
    <enumeratedValueSet variable="thy1ps-count">
      <value value="45"/>
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="thy1ns-count">
      <value value="0"/>
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-tnf-thresh">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-il1b-thresh">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="transition?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="transition-threshold">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="colvar">
      <value value="83.4"/>
      <value value="150.12"/>
      <value value="166.8"/>
      <value value="183.48"/>
      <value value="250.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="end-tick">
      <value value="365"/>
      <value value="1095"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Global Variable Sampling" repetitions="8" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>count patches with [pcolor = 128] / count patches</metric>
    <metric>count patches with [pcolor = 125] / count patches</metric>
    <metric>((count patches with [pcolor = 128] + count patches with [pcolor = 128])/ count patches) * 100</metric>
    <metric>count patches with [tissue-type = "alveolar-space" and pcolor != white] / count patches with [tissue-type = "alveolar-space"]</metric>
    <metric>count thy1ps</metric>
    <metric>count thy1ns</metric>
    <metric>fibrosis-score</metric>
    <metric>((count thy1ns / count turtles) * 100)</metric>
    <enumeratedValueSet variable="thy1ps-count">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="thy1ns-count">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-tnf-thresh">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-il1b-thresh">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="transition?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="transition-threshold">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="colvar">
      <value value="152.12"/>
      <value value="166.8"/>
      <value value="183.48"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="end-tick">
      <value value="365"/>
      <value value="1095"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prolif-rate">
      <value value="0.0162"/>
      <value value="0.018"/>
      <value value="0.0198"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="apop-rate">
      <value value="0.072"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="daily-migration-rate">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Global Variable Sampling_Limited t" repetitions="5" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count patches with [pcolor = 128] / count patches</metric>
    <metric>count patches with [pcolor = 125] / count patches</metric>
    <metric>((count patches with [pcolor = 128] + count patches with [pcolor = 128])/ count patches) * 100</metric>
    <metric>count patches with [tissue-type = "alveolar-space" and pcolor != white] / count patches with [tissue-type = "alveolar-space"]</metric>
    <metric>count thy1ps</metric>
    <metric>count thy1ns</metric>
    <metric>fibrosis-score</metric>
    <enumeratedValueSet variable="thy1ps-count">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="thy1ns-count">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-tnf-thresh">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-il1b-thresh">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="transition?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="transition-threshold">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="colvar">
      <value value="166.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="end-tick">
      <value value="365"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prolif-rate">
      <value value="0.0162"/>
      <value value="0.018"/>
      <value value="0.0198"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="apop-rate">
      <value value="0.072"/>
      <value value="0.08"/>
      <value value="0.088"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="daily-migration-rate">
      <value value="0.75"/>
      <value value="1"/>
      <value value="1.25"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Test Gen Tracking" repetitions="5" sequentialRunOrder="false" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>[(word myid) ] of turtles</metric>
    <metric>[(word psource) ] of turtles</metric>
    <metric>[(word nsource) ] of turtles</metric>
    <metric>image-title</metric>
    <enumeratedValueSet variable="thy1ps-count">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="thy1ns-count">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-tnf-thresh">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-il1b-thresh">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="transition?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="transition-threshold">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="colvar">
      <value value="166.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="end-tick">
      <value value="365"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prolif-rate">
      <value value="0.018"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="apop-rate">
      <value value="0.072"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="daily-migration-rate">
      <value value="1.25"/>
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
