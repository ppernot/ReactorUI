docker kill reactorui
docker rm reactorui

cd ~/Bureau/Titan-APSIS/
docker build -t ppernot1/reactorui:v1.5 -t ppernot1/reactorui:latest -f ReactorUI/docker/Dockerfile .

cd ~/Bureau/ChemDBTest/Projects/

docker run -d -p 3838:3838 --mount type=bind,source="$(pwd)",target=/Projects --mount type=bind,source="$(pwd)"/../ChemDBPublic,target=/ChemDBLocal --name reactorui ppernot1/reactorui
