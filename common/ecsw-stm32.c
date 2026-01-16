#include "ecsw.h"
#include "../common/stm32wrapper.h"

void setup(const enum ecsw_mode mode) {
  clock_setup(mode == ECSW_BENCHMARK ? CLOCK_BENCHMARK : CLOCK_FAST);
  gpio_setup();
  usart_setup(115200);

  SCS_DEMCR |= SCS_DEMCR_TRCENA;
  DWT_CYCCNT = 0;
  DWT_CTRL |= DWT_CTRL_CYCCNTENA;

  gpio_toggle(GPIOD, GPIO14);
}

void teardown(void) {
  gpio_toggle(GPIOD, GPIO12 | GPIO14);
  while(1);
}

void output(char *message) {
  send_USART_str((unsigned char*)message);
}

unsigned long long performance_count(void) {
  return DWT_CYCCNT;
}

char* performance_count_unit(void) {
  return "cycles";
}
