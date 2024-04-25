def tokenize_dict(token, zinc250_vocab = None):
	''' 
	Dynamically update the vocabulary or convert
	the moecules to the dictionary conversions
	'''

	# Will have to update to a class object so that 
	# entries can be added and saved
	# Built using Kaggle Zinc250K dataset
	if zinc250_vocab is None:
		zinc250_vocab = {'CC': 0,
						'C TQ ': 1,
						'CO none== 42': 2,
						'CCCCCC=': 3,
						'NC': 4,
						'CCCCO=': 5,
						'CN': 6,
						'CF': 7,
						'CNCNN=': 8,
						'CCCCCN=': 9,
						'CCCCCC': 10,
						'CN=== 43': 11,
						'CO': 12,
						'OC': 13,
						'CCCCN': 14,
						'CCCCCCN': 15,
						'CCCCCN': 16,
						'CCNCN=': 17,
						'CCl': 18,
						'CS': 19,
						'SC': 20,
						'N TQ ': 21,
						'CCCCS=': 22,
						'CBr': 23,
						'NO': 24,
						'CN none== 43': 25,
						'CN trans== 33': 26,
						'CCCCN=': 27,
						'CCCNCN=': 28,
						'CCCCC': 29,
						'CCCNN=': 30,
						'CCNCN': 31,
						'CCCCO': 32,
						'CCNCCO': 33,
						'CS none== 42': 34,
						'CCNCS=': 35,
						'NS': 36,
						'OS none== 26': 37,
						'S TQ ': 38,
						'CCNCCN': 39,
						'CCC': 40,
						'CC none== 33': 41,
						'CCCNN': 42,
						'CCCCNN=': 43,
						'SN': 44,
						'CCCC': 45,
						'NO none== 42': 46,
						'CNNCS=': 47,
						'CCCNO=': 48,
						'CNCNO=': 49,
						'CC none== 23': 50,
						'CCCCS': 51,
						'CNNCO=': 52,
						'CC none== 43': 53,
						'CCCNO': 54,
						'CCCCCO=': 55,
						'CC cis== 34': 56,
						'CO cis== 42': 57,
						'CCCCCO': 58,
						'CCCNCCN': 59,
						'CCNCO=': 60,
						'CCCCCS': 61,
						'CCNCCN=': 62,
						'CN none== 33': 63,
						'NN': 64,
						'CCCNCN': 65,
						'CCO': 66,
						'CCNCNN=': 67,
						'CO trans== 42': 68,
						'CC trans== 34': 69,
						'CCNCS': 70,
						'CNNNN=': 71,
						'CCNCCS': 72,
						'CC none== 34': 73,
						'CC trans== 33': 74,
						'CC=== 44': 75,
						'CCOCCO': 76,
						'CC none== 24': 77,
						'CCCO': 78,
						'CC=== 34': 79,
						'CCCNS': 80,
						'CCCNCO': 81,
						'CO none== 32': 82,
						'ClC': 83,
						'OS none== 24': 84,
						'CCCCCCCN': 85,
						'CCOCO': 86,
						'CS cis== 42': 87,
						'CC trans== 44': 88,
						'CCNON=': 89,
						'CCCCCCC': 90,
						'FC': 91,
						'CCNNN=': 92,
						'CNCNS=': 93,
						'CC none== 44': 94,
						'CN cis== 33': 95,
						'CN cis== 44': 96,
						'CCNNS=': 97,
						'CC cis== 44': 98,
						'CN trans== 34': 99,
						'CN none== 44': 100,
						'CCCOCCO': 101,
						'OS cis== 26': 102,
						'CN cis== 43': 103,
						'CC trans== 43': 104,
						'CCNCNN': 105,
						'CCNCO': 106,
						'CN trans== 43': 107,
						'CCNCNS': 108,
						'CS trans== 42': 109,
						'CC cis== 43': 110,
						'CCCCCCCCC': 111,
						'SO': 112,
						'CNCNCN=': 113,
						'CCNSN=': 114,
						'CC none== 32': 115,
						'CCCN': 116,
						'CCCNCCO': 117,
						'IC': 118,
						'ON': 119,
						'CCCNCCS': 120,
						'CI': 121,
						'CCCSCS': 122,
						'OS': 123,
						'CNNCS': 124,
						'CN trans== 42': 125,
						'CCCCCCO': 126,
						'CN none== 42': 127,
						'NO trans== 42': 128,
						'CNCNN': 129,
						'CCCNCS': 130,
						'BrC': 131,
						'NO cis== 42': 132,
						'CN trans== 44': 133,
						'CCNS': 134,
						'CNCNCN': 135,
						'CCCCCCS': 136,
						'CCCNS=': 137,
						'CCCOCO': 138,
						'P TQ ': 139,
						'NP none== 35': 140,
						'PN': 141,
						'PC': 142,
						'OS trans== 26': 143,
						'CCSCCS': 144,
						'CCCCNS': 145,
						'CCCCCCCC': 146,
						'OO': 147,
						'CCCCCS=': 148,
						'CCOCCS': 149,
						'CNCNCS': 150,
						'CNNCO': 151,
						'CCCCCNCS': 152,
						'CCCCNN': 153,
						'CCN': 154,
						'CCCNCO=': 155,
						'CCCC=': 156,
						'CO none== 43': 157,
						'CCCCCCNCN': 158,
						'CCCCNCCN': 159,
						'CCNNCS': 160,
						'CCCNSN': 161,
						'NN cis== 33': 162,
						'CCCCCNCN': 163,
						'CC=== 43': 164,
						'NN trans== 33': 165,
						'NN none== 24': 166,
						'NN none== 43': 167,
						'CC none== 42': 168,
						'NP': 169,
						'OP none== 25': 170,
						'CCCNNN=': 171,
						'CCCCCCCCN': 172,
						'CCCCCNCCN': 173,
						'CCNNO=': 174,
						'CCCCCCCNN': 175,
						'CCCSS': 176,
						'CNCNO': 177,
						'CN none== 34': 178,
						'CCCCNCN': 179,
						'CCOCCOCCOCCOCCO': 180,
						'NS trans== 34': 181,
						'CCCNCNCNN': 182,
						'CCNCNCNCN': 183,
						'CCCNCNNCN': 184,
						'CCCOCCS': 185,
						'NN none== 34': 186,
						'NN none== 42': 187,
						'CCCCNCO': 188,
						'CCCCNCNCO': 189,
						'CC cis== 33': 190,
						'CCCCCCNCCN': 191,
						'NN=== 43': 192,
						'OP': 193,
						'CCNNN': 194,
						'CN cis== 42': 195,
						'CNCNS': 196,
						'CCOCS': 197,
						'CP': 198,
						'PO': 199,
						'NO none== 32': 200,
						'SS': 201,
						'CCCS': 202,
						'CN cis== 34': 203,
						'CCNSS=': 204,
						'CCCCCCC=': 205,
						'NS none== 36': 206,
						'CCCNCNCCO': 207,
						'CCCCNCCCN': 208,
						'CCCSCCS': 209,
						'CCCCCCCCCCC': 210,
						'CCCCOCCO': 211,
						'CCSCS': 212,
						'CCOCCOCCOCCO': 213,
						'CCCCCNN': 214,
						'NN=== 34': 215,
						'CCNNCP': 216,
						'CCCCCCCCCCCC': 217,
						'CCCCCCCCNCCN': 218,
						'CCCCCCCNCCCN': 219,
						'CCCNNS': 220,
						'CCNCCS=': 221,
						'CCCCCCCO': 222,
						'CCOCCOCCOCCOCCOCCOCCOCCO': 223,
						'NN cis== 34': 224,
						'CCNCCO=': 225,
						'CCCCCCCCO': 226,
						'CCCCCOCO': 227,
						'CCCSS=': 228,
						'CCCNNCO': 229,
						'CCCCNPN': 230,
						'CNNNN': 231,
						'CCCCCCCCCCCCCO': 232,
						'CCCCNCN=': 233,
						'CCCCC=': 234,
						'CO cis== 32': 235,
						'CCCCCCCCCC': 236,
						'CCCNCS=': 237,
						'CCCNCNN': 238,
						'CNN': 239,
						'CNCN': 240,
						'CCCCNO': 241,
						'CCCNCCCN': 242,
						'NN none== 33': 243,
						'CNCSS=': 244,
						'CCCNCCCO': 245,
						'CSCSCS': 246,
						'CO trans== 32': 247,
						'CCNNP': 248,
						'PS none== 52': 249,
						'OS trans== 24': 250,
						'CCCCCCN=': 251,
						'CCCCNCNCN': 252,
						'CCCCCNCNN': 253,
						'CCOCS=': 254,
						'OP trans== 25': 255,
						'CCCOPO': 256,
						'OS cis== 24': 257,
						'CCNCCNCS': 258,
						'CCCCNCCNCN': 259,
						'CCCCCNCCNN': 260,
						'CCCNCCCNO': 261,
						'CS none== 32': 262,
						'CCCCNCS': 263,
						'CCCCNCCO': 264,
						'CCCCCNCCO': 265,
						'CCCCCCCCNN': 266,
						'NN none== 44': 267,
						'CCCCOCCOCCNCCOCCOCCNCCOCCO': 268,
						'NS none== 24': 269,
						'CCCCCCCCCCCCCCCC': 270,
						'CCCCCCCCCCCCCC': 271,
						'CCCCNCCCO': 272,
						'CCCNP': 273,
						'CS none== 43': 274,
						'CCNNS': 275,
						'CCNCCNN': 276,
						'CCCOS': 277,
						'CCCCCCCCCCN': 278,
						'SCl': 279,
						'CCOCCOCCOCCOCCOCCS': 280,
						'CCCCOCCCO': 281,
						'CCNCCOCCOCCO': 282,
						'CCCCCNCNO': 283,
						'CCS': 284,
						'CNNCNN=': 285,
						'CCCCCCNO': 286,
						'CNCS': 287,
						'CCCCCCNN': 288,
						'CO trans== 43': 289,
						'CCSCS=': 290,
						'CCCNNCS': 291,
						'CNNCNN': 292,
						'SP': 293,
						'CCSCCS=': 294,
						'CCNNCNN': 295,
						'NN trans== 34': 296,
						'NN trans== 42': 297,
						'CCCOS=': 298,
						'CNNNS=': 299,
						'SS none== 62': 300,
						'CNCNCO': 301,
						'CS cis== 43': 302,
						'CCNOCO': 303,
						'CNCNCO=': 304,
						'CCNCCNS': 305,
						'CCSSCCSS': 306,
						'CCCCOCCCCO': 307,
						'CCCCCOCCO': 308,
						'CCCNSN=': 309,
						'CCCNCCOCN': 310,
						'CCNSN': 311,
						'CCNSCO': 312,
						'CCSCCSCCS': 313,
						'FN': 314,
						'CCCNPN': 315,
						'CCNPO': 316,
						'CCCNCCNCO': 317,
						'NCl': 318,
						'NBr': 319,
						'CCSNS': 320,
						'CCOSO': 321,
						'CCCCCCCCSCCCNCS': 322,
						'CCCCCCNCSCCCCCCS': 323,
						'CCNNCCSS': 324,
						'NNNN': 325,
						'CCOPO': 326,
						'CCCCCCCOCCCCNCCCO': 327,
						'CCCCCCCNO': 328,
						'CCCOCSCCS': 329,
						'CCNNCS=': 330,
						'CN none== 32': 331,
						'CCCCCNCCOCCNCCOCCN': 332,
						'CCCSCCSCCCSCCS': 333,
						'CCCCCNCON': 334,
						'CCCNPO': 335,
						'CNCSS': 336,
						'CCCCOCO': 337,
						'CCCCCCCCS': 338,
						'CCNNCO': 339,
						'FO': 340,
						'CCCCCCCS': 341,
						'CCCCCCCCCN': 342,
						'CNOCNO': 343,
						'CCCCNCCCCN': 344,
						'CCNCCOCCOCCNCCOCCO': 345,
						'CNCOCS': 346,
						'CCNCNO': 347,
						'CCNCCOCCOCCOCCOCCO': 348,
						'CCNSCCSO': 349,
						'CCCNCCN=': 350,
						'CCNCP': 351,
						'CCCCCCCCCCCO': 352,
						'CCCCCNCO': 353,
						'CCNCNP': 354,
						'NF': 355,
						'CCOCCOCCOCCOCCOCCO': 356,
						'NO trans== 32': 357,
						'NS none== 26': 358,
						'CCNNCO=': 359,
						'PS': 360,
						'CCNCCNCCNCCNCCNCCN': 361,
						'CCCCNP': 362,
						'CNCNPN': 363,
						'CCNCCOCCOCCOCCO': 364,
						'NN none== 23': 365,
						'CCOCPCO': 366,
						'CNOSN': 367,
						'CCCCCOCOO': 368,
						'CCNCOCCS': 369}
		

	token_count = len(list(zinc250_vocab.keys()))

	# if token in zinc250_vocab:
	# 	conversion = zinc250_vocab[token]

	# if token not in zinc250_vocab:
	# 	print()
		# token_count          += 1
		# zinc250_vocab[token] = token_count 
		# conversion           = token_count

	return zinc250_vocab