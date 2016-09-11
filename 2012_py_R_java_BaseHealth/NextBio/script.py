#! /usr/bin/env python
import json,pymongo
 
# mongo_cnxn_string = "mongodb://adm1n:adm1npw0rd@54.244.24.168:27017/admin"
# conn = pymongo.Connection(mongo_cnxn_string, 47954)
# db = conn["crm"]
# 
# mongo_cnxn_string2 = "mongodb://adm1n:adm1npw0rd@54.244.24.253:27017/admin"
# conn2 = pymongo.Connection(mongo_cnxn_string2, 27017)
# db2 = conn2["genophen30"]

######### tiger
mongo_cnxn_string2 = "mongodb://adm1n:adm1npw0rd@54.245.239.40:27017/admin"
conn2 = pymongo.Connection(mongo_cnxn_string2, 27017)
db2 = conn2["genophen30"]


def findSashaAuthentication():
    print "findSashaAuthentication"
    members = db2["sasha.authentication"].find()
    pidToDelete = []  
    for m in members:
        pid = m['static']['gid']
        f = False
        person = db['persons'].find({}, { "info.pid" : 1}) 
        for p in person:
#             print p['info']['pid'], pid
            if pid == p['info']['pid']:
#                 print "found" , p
                f = True
                break
        if f == False:
            if pid != "0" and pid !="00acc4a9-489d-4d62-8630-766c8611160f":
                print "not found", pid
#                 db2["sasha.authentication"].remove({"static.gid":pid})

def deleteMembers( pidsToDelete ):
    print "Deleting members"
    for pid in pidsToDelete:
        print "Delete " + pid
        #db["members"].remove( { "pid" : pid } )

def findAllMembers():
    print "Finding members"
    members = db["members"].find()
    return members

def findAllPhysicians():
    print "Finding findAllPhysicians"
    members = db["physicians"].find()
    return members

def findPersonsForMembers( members ):
    print "Finding findPersonsForMembers"
    deleteList = []
    for member in members:
#         print member["pid"]
        if db["persons"].find( { "info.pid" : member["pid"] } ).count() == 0 :
            deleteList.append( member["pid"])
    return deleteList


def findAllCaregroup(  ):
    print "Finding caregroups"
    careGroups  = db["careGroups"].find()
    return careGroups 

def findPersonsForCaregroup( members ):
    print "Finding findPersonsForCaregroup"
    deleteByMember = []
    deleteByPhysician = []
    for member in members:
#         print  member["memberPid"],  member["careGroupId"]
        if db["persons"].find( { "info.pid" : member["memberPid"] } ).count() == 0 :
            print "person delete" , member["memberPid"]
            deleteByMember.append( member["memberPid"])
        if db["members"].find( { "pid" : member["memberPid"] } ).count() == 0 :
            print "members delete", member["memberPid"]
            deleteByMember.append( member["memberPid"])
        if db["persons"].find( { "info.pid" : member["careGroupId"] } ).count() == 0 :
            print "physisican not in persondelete", member["careGroupId"]
            deleteByPhysician.append( member["careGroupId"])
        if db["physicians"].find( { "pid" : member["careGroupId"] } ).count() == 0 :
            print "physisican delete", member["careGroupId"]
            deleteByPhysician.append( member["careGroupId"])
        
    return deleteByMember, deleteByPhysician

def deleteCargroups( pidMember, pidPhysician  ):
    print "deleteCargroups"
    for pid in pidMember:
        print "Delete member" + pid
        db["careGroups"].remove( { "memberPid" : pid } )
    for pid in pidPhysician:
        print "Delete physician" + pid
        db["careGroups"].remove( { "careGroupId" : pid } )

def findAllEntrepriseMapping(  ):
    print "Finding enterpriseGroupMapping"
    enterpriseGroupMapping  = db["enterpriseGroupMapping"].find()
    return enterpriseGroupMapping 

def findPhysicianForEntrepriseMapping( entreprises ):
    print "Finding findPersonsForCaregroup"
    deleteByPhysician = []
    for entreprise in entreprises:
        if db["persons"].find( { "info.pid" : entreprise["physicianPid"] } ).count() == 0 :
            deleteByPhysician.append( entreprise["physicianPid"])
        if db["physicians"].find( { "pid" : entreprise["physicianPid"] } ).count() == 0 :
            deleteByPhysician.append( entreprise["physicianPid"])

def findAllEntrepriseAdmin(  ):
    print "Finding findAllEntrepriseAdmin"
    enterpriseGroupAdmins  = db["enterpriseGroupAdmins"].find()
    return enterpriseGroupAdmins 

def findPersonForEntrepriseAdmin( entreprises ):
    print "Finding findPersonsForCaregroup"
    deleteByAdmin = []
    for entreprise in entreprises:
        if db["persons"].find( { "info.pid" : entreprise["pid"] } ).count() == 0 :
            deleteByAdmin.append( entreprise["pid"])

def getLostSnpInGeneticsRep():
    
    snp = conn2["nextbio"]['dbsnp'].find()
    dbsnp = []
    for doc in snp:
        if doc['dbsnp'] not in dbsnp:
            dbsnp.append(doc['dbsnp'])
    genetics_data  = db2["genetics.rep"].find()
    new_snp = []
    for doc in genetics_data:
        for data in doc['snps']:
            if data['snp'] not in dbsnp and data['snp'] not in new_snp:
                new_snp.append(data['snp'])
    print json.dumps(new_snp)
         

getLostSnpInGeneticsRep()

# findSashaAuthentication()  

# physicians = findAllPhysicians()
# pidsToDelete = findPersonsForMembers( physicians )
# print pidsToDelete

# members = findAllMembers()
# pidsToDelete = findPersonsForMembers( members )
# print pidsToDelete
# deleteMembers( pidsToDelete )

# careGroups=findAllCaregroup()
# pidMember, pidPhysician = findPersonsForCaregroup (careGroups)
# print pidMember, pidPhysician
# # deleteCargroups( pidMember, pidPhysician )
# 
# entrepriseMapping = findAllEntrepriseMapping()
# pidPhysician = findPhysicianForEntrepriseMapping(entrepriseMapping)
# print pidPhysician

# entrepriseAdmin = findAllEntrepriseAdmin()
# pidAdmin= findPersonForEntrepriseAdmin(entrepriseAdmin)
# print pidAdmin
