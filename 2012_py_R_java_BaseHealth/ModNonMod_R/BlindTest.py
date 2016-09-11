# in order to make blind testing easier for Algorithm Group and Doctors
import MySQLdb, random

class BlindTest():
     def __init__(self, id):
         self.Did     = id #DiseaseID
         self.factors = {}
         self.units   = {}
         
     
     def DBconnect(self):
         db = MySQLdb.connect("192.168.2.104","journal","journal","AlgoPhen" )
         cursor = db.cursor()
         sql = '''select Q.qtype, I.IDName, I.unit, Q.Gender, M.isFactor, M.Categories, M.MinVal, M.MaxVal 
         from ModNmod M, Qmap Q, IDmap I 
         where I.IsRiskFactorID='yes' and Q.RiskFactorID=I.IDNo and 
         M.RiskID=I.IDNo and M.Genophen='yes' and 
         M.DiseaseID=%s and Q.DiseaseID like '%%%s%%' order by Q.QID''' % (self.Did, self.Did)
         
         cursor.execute(sql)
         results = cursor.fetchall()
         
         for row in results:
             self.units[row[1]] = row[2]
             if row[0] in self.factors.keys():
                 if row[4]=='yes':
                     self.factors[row[0]][row[1]] = row[5].split(',')
                 else:
                     self.factors[row[0]][row[1]] = range(row[6], row[7]) 
             else:
                 if row[4]=='yes':
                     self.factors[row[0]] = {row[1]:row[5].split(',')}
                 else:
                     self.factors[row[0]] = {row[1]:range(row[6], row[7]) }
                 
         
         for k, v in self.factors.items():
             print "<< ", k, " >>"
             for kk, vv in v.items():
                 print "\t< ", kk, " >"
                 print "\t\t", vv
                 print "\t\t", self.units[kk]
                 print ""
                  
    
    
     def randomize(self):
         # randomly created profiles
         random.seed()
         
         bible = {1:'T2D', 4:'CRC', 5:'Alz', 6:'HD', 7:'Strk', 8:'DP', 9:'Hyp', 101:'Anx' , 102:'Bcp' , 105:'Hrl' ,113:'Edy'}
         file = open("/home/pouria3/Dropbox/EngExl/%s/%s_BlindTest.csv"%(bible[self.Did], bible[self.Did]), 'w')
         profile = ''
         
         for k in self.units.keys():
             profile += k
             if self.units[k] is not None:
                 profile += ' (' + self.units[k] + ')'  
             for num in range(0, 15):
                 for k1,v1 in self.factors.items():
                     if k in v1.keys():
                         vv = v1[k]
                         profile += ', ' + str(vv[random.randint(0, len(vv)-1)]) 
             profile += '\n'        
         print profile
         file.write(profile)
                 
#         for k, v in self.factors.items():
#             profile += '<dive align="left">\n'
#             profile += "<h1> " + k + " </h1>" + "\n"
#             for kk, vv in v.items():
                 #profile += '<dive align="center">\n'
#                 profile += "<h2>&nbsp; " + kk + "</h2>" + "</div>\n"
#                 profile += "<p>&nbsp;&nbsp" + str(vv[random.randint(0, len(vv)-1)]) + " : " + str(self.units[kk]) + "</p>" + "\n"
         
#         file.write(profile)
                 

        
         
         
         
         
         
         
         
         
         
if __name__ == "__main__":
    b = BlindTest(9)
    b.DBconnect()
    b.randomize()
    
    
         
             

         
    