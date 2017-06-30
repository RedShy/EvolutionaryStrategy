This repository contains some implementations of Evolution Strategy for computing the MPED for two given strings.

Usage is "mped $heuristic π_1 π_2" where $heuristic is the specific parameter for executing a particular heuristic. A list of parameters and their corresponding heuristic is listed below.

| First Header  | Second Header |
| ------------- | ------------- |
| Content Cell  | Content Cell  |
| Content Cell  | Content Cell  |

(1+1)-ES 						           	"es-one-one"
(1+1)-ES Simple Random Restart 	"es-one-one-srs"
(μ+1)-ES 					            	"es-wp"
(μ+λ)-ES 						            "es"
(μ,λ)-ES 						            "es-comma"
Bruteforce 						          "ex"
Fast Bruteforce 				        "exf"
Random Search 					        "random"
Hill climbing 					        "hc"
