#!/bin/bash
#Program:
#	use 'read' command to get users' name and print out
#History:
#2012/10/02 YuanZheCSYZ First release

read -p "Please input your first name:" firstname
read -p "Please input your last name:" lastname
echo -e "\nYou full name is: $firstname $lastname"
