 #include <ros.h>
#include <std_msgs/Int8MultiArray.h>

const int thr1 = A2; //assign correct pin
const int thr2 = A3; //assign correct pin
const int thr3 = A5; //assign correct pin
const int thr4 = A4; //assign correct pin
const int thr5 = A0; //assign correct pin
const int thr6 = A1; //assign correct pin
const int thr7 = 7; //assign correct pin
const int thr8 = 6; //assign correct pin

ros::NodeHandle  nh;

void thrusters_callback(const std_msgs::Int8MultiArray& arr){
  for(int i = 0; i < 8; ++i){
    switch(i){
      case 0:
        if(arr.data[i] > 0){
          digitalWrite(thr1, HIGH);
        }else{
          digitalWrite(thr1, LOW);
        }
        break;
      case 1:
        if(arr.data[i] > 0){
          digitalWrite(thr2, HIGH);
        }else{
          digitalWrite(thr2, LOW);
        }
        break;
      case 2:
        if(arr.data[i] > 0){
          digitalWrite(thr3, HIGH);
        }else{
          digitalWrite(thr3, LOW);
        }
        break;
      case 3:
        if(arr.data[i] > 0){
          digitalWrite(thr4, HIGH);
        }else{
          digitalWrite(thr4, LOW);
        }
        break;
      case 4:
        if(arr.data[i] > 0){
          digitalWrite(thr5, HIGH);
        }else{
          digitalWrite(thr5, LOW);
        }
        break;
      case 5:
        if(arr.data[i] > 0){
          digitalWrite(thr6, HIGH);
        }else{
          digitalWrite(thr6, LOW);
        }
        break;
      case 6:
        if(arr.data[i] > 0){
          digitalWrite(thr7, HIGH);
        }else{
          digitalWrite(thr7, LOW);
        }
        break;
      case 7:
        if(arr.data[i] > 0){
          digitalWrite(thr8, HIGH);
        }else{
          digitalWrite(thr8, LOW);
        }
        break;
    }
  }

}

ros::Subscriber<std_msgs::Int8MultiArray> sub("/thrusters", &thrusters_callback );

void setup() {
  pinMode(thr1, OUTPUT);
  pinMode(thr2, OUTPUT);
  pinMode(thr3, OUTPUT);
  pinMode(thr4, OUTPUT);
  pinMode(thr5, OUTPUT);
  pinMode(thr6, OUTPUT);
  pinMode(thr7, OUTPUT);
  pinMode(thr8, OUTPUT);
  nh.initNode();
  nh.subscribe(sub);

}

void loop() {
  nh.spinOnce();
}
