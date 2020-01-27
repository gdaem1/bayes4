library(shiny)
library(plotly)

shinyUI(
  fluidPage(
    titlePanel('Bayesian A/B/n test for ARPU'),
    sidebarLayout(
      sidebarPanel(
        numericInput(
          'rev_A',
          'Total Revenue for Control:',
          min = 0,
          max = 1e6,
          value = 4675.13
        ),
        numericInput(
          'success_A',
          'Control Conversions:',
          min = 0,
          max = 1e6,
          value = 405
        ),
        numericInput(
          'total_A',
          'Total Control Players:',
          min = 0,
          max = 1e6,
          value = 32597
        ),
        numericInput(
          'retention_A',
          'retention A:',
          min = 0,
          max = 1e6,
          value = 200
        ),
        numericInput(
          'rev_B',
          'Total Revenue Treatment 1:',
          min = 0,
          max = 1e6,
          value = 2534.80
        ),
        numericInput(
          'success_B',
          'Converted Players Treatment 1:',
          min = 0,
          max = 1e6,
          value = 212
        ),
        numericInput(
          'total_B',
          'Total Players Treatment 1:',
          min = 0,
          max = 1e6,
          value = 16437
        ),
        numericInput(
          'retention_B',
          'retention B:',
          min = 0,
          max = 1e6,
          value = 230
        ),
        numericInput(
          'rev_C',
          'Total Revenue Condition 2:',
          min = 0,
          max = 1e6,
          value = 2606.49
        ),
        numericInput(
          'success_C',
          'Converted Players Condition 2:',
          min = 0,
          max = 1e6,
          value = 227
        ),
        numericInput(
          'total_C',
          'Total Players in Condition 2:',
          min = 0,
          max = 1e6,
          value = 16199
        ),
        numericInput(
          'retention_C',
          'retention C:',
          min = 0,
          max = 1e6,
          value = 230
        ),
        numericInput(
          'rev_D',
          'Total Revenue Condition 3:',
          min = 0,
          max = 1e6,
          value = 2606.49
        ),
        numericInput(
          'success_D',
          'Converted Players Condition 3:',
          min = 0,
          max = 1e6,
          value = 227
        ),
        numericInput(
          'total_D',
          'Total Players in Condition 3:',
          min = 0,
          max = 1e6,
          value = 16199
        ),
        numericInput(
          'retention_D',
          'retention D:',
          min = 0,
          max = 1e6,
          value = 230
        ),
        numericInput(
          'sim_sample',
          'Monte-Carlo Simulation Samples:',
          min = 2,
          max = 1e7,
          value = 1e5
        ),
        actionButton(
          'button',
          'Calculate'
        ),
        hr(),
        tags$div(
          class='header', checked = NA,
          tags$p(
            'Bayesian test for A/B/n ARPU calculator. Modified for internal use'
          ),
          tags$br(),
          tags$a(
            href = 'https://github.com/Vidogreg/bayes-ab-testing/tree/master/bayes-arpu-test',
            'Original Code by V. Gregor'
          ),
          tags$p('')
        )
      ),
      mainPanel(
        tableOutput('table1'),
        tableOutput('table2'),
        tableOutput('table3'),
        tableOutput('table4')
      )
    )
  )
)
