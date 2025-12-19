const int thr1 = A2;
const int thr2 = A3;
const int thr3 = A5;
const int thr4 = A4;
const int thr5 = A0;
const int thr6 = A1;
const int thr7 = 7;
const int thr8 = 6;

const int NUM_THRUSTERS = 8;

int thrusterPins[NUM_THRUSTERS] = {
  thr1, thr2, thr3, thr4, thr5, thr6, thr7, thr8
};

void setup() {
  Serial.begin(9600);  // Make sure your PC matches this baud rate

  for (int i = 0; i < NUM_THRUSTERS; i++) {
    pinMode(thrusterPins[i], OUTPUT);
    digitalWrite(thrusterPins[i], LOW);  // Ensure all thrusters start off
  }

  Serial.println("Thruster board started...");
}

void loop() {
  // Wait until 8 bytes are available
  if (Serial.available() >= NUM_THRUSTERS) {
    int8_t values[NUM_THRUSTERS];

    for (int i = 0; i < NUM_THRUSTERS; i++) {
      int raw = Serial.read();  // Read one byte
      values[i] = (int8_t)raw;  // Convert to signed int8
    }

    char buffer[128];
    sprintf(buffer, "New data available: %d %d %d %d %d %d %d %d", values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7]);
    Serial.println(buffer);

    // Process and set each thruster pin
    for (int i = 0; i < NUM_THRUSTERS; i++) {
      if (values[i] > 0) {
        digitalWrite(thrusterPins[i], HIGH);
      } else {
        digitalWrite(thrusterPins[i], LOW);
      }
    }
  }
}
